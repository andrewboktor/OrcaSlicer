#include "SpiralVase.hpp"
#include "GCode.hpp"
#include <sstream>
#include <cmath>
#include <limits>

namespace Slic3r {

float distance(SpiralPoint a, SpiralPoint b) {
    return sqrt(pow(a.x-b.x,2)+pow(a.y-b.y, 2));
}

/** Find nearest point to p from population, returns the index*/
long nearest(SpiralPoint p, std::vector<SpiralPoint> *population)
{
    float min     = std::numeric_limits<float>::max();
    long  min_idx = -1;
    for (unsigned long i = 0; i < population->size(); i++) {
        float dist = distance(population->at(i), p);
        if (dist < min) {
            min     = dist;
            min_idx = i;
        }
    }
    return min_idx;
}

SpiralPoint subtract(SpiralPoint a, SpiralPoint b) {
    return SpiralPoint(a.x-b.x, a.y-b.y);
}

SpiralPoint add(SpiralPoint a, SpiralPoint b) {
    return SpiralPoint(a.x+b.x, a.y+b.y);
}

SpiralPoint scale(SpiralPoint a, float factor){
    return SpiralPoint(a.x*factor, a.y*factor);
}

float dot(SpiralPoint a, SpiralPoint b) {
    return a.x*b.x+a.y*b.y;
}

SpiralPoint nearest_point_on_line(SpiralPoint c, SpiralPoint a, SpiralPoint b, float& dist) {
    SpiralPoint ab = subtract(b, a);
    SpiralPoint ca = subtract(c, a);
    float t = dot(ca, ab)/dot(ab,ab);
    t=t>1?1:t;
    t=t<0?0:t;
    SpiralPoint closest= SpiralPoint(add(a, scale(ab, t)));
    dist = distance(c, closest);
    return closest;
}

SpiralPoint nearest_point_on_polygon(SpiralPoint p, std::vector<SpiralPoint>* polygon, bool& found, float& dist) {
    if(polygon->size()<2) {
        found=false;
        return SpiralPoint(0,0);
    }
    float min     = std::numeric_limits<float>::max();
    SpiralPoint closest(0,0);
    for(unsigned long i=0; i<polygon->size()-1; i++) {
        float currentDist=0;
        SpiralPoint current = nearest_point_on_line(p, polygon->at(i), polygon->at(i+1), currentDist);
        if(currentDist<min) {
            min=currentDist;
            closest=current;
            found=true;
        }
    }
    dist=min;
    return closest;

}

std::string SpiralVase::process_layer(const std::string &gcode, bool last_layer)
{
    /*  This post-processor relies on several assumptions:
        - all layers are processed through it, including those that are not supposed
          to be transformed, in order to update the reader with the XY positions
        - each call to this method includes a full layer, with a single Z move
          at the beginning
        - each layer is composed by suitable geometry (i.e. a single complete loop)
        - loops were not clipped before calling this method  */
    
    // If we're not going to modify G-code, just feed it to the reader
    // in order to update positions.
    if (! m_enabled) {
        m_reader.parse_buffer(gcode);
        return gcode;
    }
    
    // Get total XY length for this layer by summing all extrusion moves.
    float total_layer_length = 0;
    float layer_height = 0;
    float z = 0.f;
    
    {
        //FIXME Performance warning: This copies the GCodeConfig of the reader.
        GCodeReader r = m_reader;  // clone
        bool set_z = false;
        r.parse_buffer(gcode, [&total_layer_length, &layer_height, &z, &set_z]
            (GCodeReader &reader, const GCodeReader::GCodeLine &line) {
            if (line.cmd_is("G1")) {
                if (line.extruding(reader)) {
                    total_layer_length += line.dist_XY(reader);
                } else if (line.has(Z)) {
                    layer_height += line.dist_Z(reader);
                    if (!set_z) {
                        z = line.new_Z(reader);
                        set_z = true;
                    }
                }
            }
        });
    }

    // Remove layer height from initial Z.
    z -= layer_height;

    std::vector<SpiralPoint>* current_layer = new std::vector<SpiralPoint>();
    std::vector<SpiralPoint>* previous_layer = m_previous_layer;

    bool smooth_spiral = m_smooth_spiral;
    std::string new_gcode;
    std::string transition_gcode;
    float max_xy_smoothing = 1; // Made up threshold to prevent snapping to points too far away, Cura uses (2*line_width)^2 but I can't figure out how to get that config for this layer
    //FIXME Tapering of the transition layer only works reliably with relative extruder distances.
    // For absolute extruder distances it will be switched off.
    // Tapering the absolute extruder distances requires to process every extrusion value after the first transition
    // layer.
    bool  transition = m_transition_layer && m_config.use_relative_e_distances.value;
    float len = 0.f;
    SpiralPoint last_point = previous_layer != NULL && previous_layer->size() >0? previous_layer->at(previous_layer->size()-1): SpiralPoint(0,0);
    m_reader.parse_buffer(gcode, [&new_gcode, &z, total_layer_length, layer_height, transition, &len, &current_layer, &previous_layer, &transition_gcode, &last_layer, &smooth_spiral, &max_xy_smoothing, &last_point]
        (GCodeReader &reader, GCodeReader::GCodeLine line) {
        if (line.cmd_is("G1")) {
            if (line.has_z()) {
                // If this is the initial Z move of the layer, replace it with a
                // (redundant) move to the last Z of previous layer.
                line.set(reader, Z, z);
                new_gcode += line.raw() + '\n';
                return;
            } else {
                float dist_XY = line.dist_XY(reader);
                if (dist_XY > 0) {
                    // horizontal move
                    if (line.has_e()) {
                        len += dist_XY;
                        float factor = len / total_layer_length;
                        if (transition)
                            // Transition layer, interpolate the amount of extrusion from zero to the final value.
                            line.set(reader, E, line.e() * factor);
                        else if (last_layer) {
                            // We want the last layer to ramp down extrusion, but without changing z height!
                            // So clone the line before we mess with its Z and duplicate it into a new layer that ramps down E
                            // We add this new layer at the very end
                            GCodeReader::GCodeLine transitionLine(line);
                            transitionLine.set(reader, E, line.value(E) * (1 - factor));
                            transition_gcode += transitionLine.raw() + '\n';
                        }
                        // This line is the core of Spiral Vase mode, ramp up the Z smoothly
                        line.set(reader, Z, z + factor * layer_height);
                        if (smooth_spiral) {
                            // Now we also need to try to interpolate X and Y
                            SpiralPoint p(line.x(), line.y()); // Get current x/y coordinates
                            current_layer->push_back(p);       // Store that point for later use on the next layer

                            if (previous_layer != NULL) {
                                bool        found    = false;
                                float       dist     = 0;
                                SpiralPoint nearestp = nearest_point_on_polygon(p, previous_layer, found, dist);
                                if (found && dist < max_xy_smoothing) {
                                    // Interpolate between the point on this layer and the point on the previous layer
                                    SpiralPoint target = add(scale(nearestp, 1 - factor), scale(p, factor));
                                    line.set(reader, X, target.x);
                                    line.set(reader, Y, target.y);
                                    // We need to figure out the distance of this new line!
                                    float modified_dist_XY = distance(last_point, target);
                                    line.set(reader, E,
                                             line.e() * modified_dist_XY / dist_XY); // Scale the extrusion amount according to change in length
                                    last_point = target;
                                } else {
                                    last_point = p;
                                }
                            }
                        }
                        new_gcode += line.raw() + '\n';
                    }
                    return;
                    /*  Skip travel moves: the move to first perimeter point will
                        cause a visible seam when loops are not aligned in XY; by skipping
                        it we blend the first loop move in the XY plane (although the smoothness
                        of such blend depend on how long the first segment is; maybe we should
                        enforce some minimum length?).  */
                }
            }
        }
        new_gcode += line.raw() + '\n';
        if(last_layer) {
            transition_gcode += line.raw() + '\n';
        }
    });

    delete m_previous_layer;
    m_previous_layer = current_layer;
    
    return new_gcode + transition_gcode;
}

}
