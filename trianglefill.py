# (C) 2015 Petteri Aimonen <jpa (a) git.mail.kapsi.fi>
#
# ***** BEGIN GPL LICENSE BLOCK *****
#
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation,
# Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301, USA.
#
# ***** END GPL LICENCE BLOCK *****

bl_info = {
    "name": "Triangle fill",
    "author": "Petteri Aimonen",
    "version": (1, 0),
    "blender": (2, 76, 0),
    "location": "",
    "description": "Fills the selected face with evenly sized triangles",
    "warning": "",
    "wiki_url": "",
    "category": "",
    }

import bpy
import bmesh
import itertools
from mathutils import Vector, geometry
import math
import time

class Plane:
    '''Represents a plane in 3 dimensional space and provides conversion
    to 2d coordinates on the plane and back to 3d.
    '''
    def __init__(self, v1, v2, v3):
        self.orig = v1
        self.x_dir = (v2 - self.orig).normalized()
        self.y_dir = v3 - self.orig
        self.y_dir = (self.y_dir - self.y_dir.project(self.x_dir)).normalized()
    
    def to_2d(self, p):
        p = p - self.orig
        x = p.dot(self.x_dir)
        y = p.dot(self.y_dir)
        return Vector([x,y])
    
    def to_3d(self, p):
        return self.orig + self.x_dir * p.x + self.y_dir * p.y

def iter_tuples(lst, n):
    '''Iterate list as n-sized tuples.
    For example iter_tuples((1,2,3), 2) => ((1,2), (2,3), (3,1))
    '''
    lst2 = lst + lst[:n]
    for i in range(len(lst)):
        yield lst2[i:i+n]

def distance_point_segment(point, v1, v2):
    '''Compute distance of a point from a line segment.'''
    x, d = geometry.intersect_point_line(point, v1, v2)
    if d <= 0:
        return v1, (point - v1).magnitude
    elif d >= 1.0:
        return v2, (point - v2).magnitude
    else:
        return x, (point - x).magnitude

class Polygon:
    '''Represents a polygon that lies in a flat plane.'''
    def __init__(self, points):
        # Find a set of three vertices that are not on same line
        for v1, v2, v3 in iter_tuples(points, 3):
            d = abs((v1 - v2).normalized().dot((v2 - v3).normalized()))
            if d < 0.99:
                self.plane = Plane(v1, v2, v3)
                break
        else:
            raise Exception("Could not find non-colinear pair of edges")
        
        self.points = [self.plane.to_2d(p) for p in points]
        self.min_x = min([p.x for p in self.points])
        self.min_y = min([p.y for p in self.points])
        self.max_x = max([p.x for p in self.points])
        self.max_y = max([p.y for p in self.points])
    
    def shortest_edge_length(self):
        '''Returns length of the shortest polygon edge'''
        return min((p1 - p2).magnitude for p1, p2 in iter_tuples(self.points, 2))
    
    def average_edge_length(self):
        '''Returns average length of polygon edges'''
        return sum((p1 - p2).magnitude for p1, p2 in iter_tuples(self.points, 2)) / len(self.points)
    
    def contains(self, point):
        '''Returns true if the 2D point is inside the polygon.'''
        # From: http://blenderartists.org/forum/showthread.php?229303-Point-in-Polygon-Script
        x, y = point[:2]
        n = len(self.points)
        inside = False
        p1x, p1y = self.points[0][:2]
        for i in range(n + 1):
            p2x, p2y = self.points[i % n][:2]
            if y > min(p1y,p2y) and y <= max(p1y,p2y) and x <= max(p1x,p2x):
                if p1y != p2y:
                    xinters = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                if p1x == p2x or x <= xinters:
                    inside = not inside
            p1x,p1y = p2x,p2y
        return inside
    
    def edge_distance(self, point):
        '''Given a 2D point, returns the closest point on polygon edge and the distance to it.'''
        best = None
        for v1, v2 in iter_tuples(self.points, 2):
            v, d = distance_point_segment(point, v1, v2)
            if best is None or d < best[1]:
                best = v, d
        return best

def sortface(pointidx, face):
    '''Sort the points in a face so that duplicates can be eliminated.'''
    face = tuple(sorted(face))
    p1, p2, p3 = [pointidx[p] for p in face]
    if (p2 - p1).cross(p3 - p2) < 0:
        face = (face[0], face[2], face[1])
    return face

def sortedge(edge):
    return tuple(sorted(edge))

class TriangleFill(bpy.types.Operator):
    bl_idname = "mesh.triangle_fill"
    bl_label = "Triangle Fill"
    bl_options = {'REGISTER', 'UNDO'}

    edge_length = bpy.props.FloatProperty(
        name = "Edge length",
        default = 100.0,
        min = 0.1,
        max = 500.0,
        subtype='PERCENTAGE',
        unit='LENGTH',
        description='Adjust target edge length as percentage of current average edge length.'
    )
    
    only_edges = bpy.props.BoolProperty(
        name = "Add only edges",
        default = False
    )

    def execute(self, context):
        start = time.time()
    
        # Load the bmesh from current edit mode object
        obj = bpy.context.object
        bm = bmesh.from_edit_mesh(obj.data)
        face = bm.faces.active # TODO: Allow selecting multiple faces
        
        # Compute target length of each edge
        poly = Polygon([v.co for v in face.verts])
        targetlen = self.edge_length * poly.average_edge_length() / 100.0
        
        # Check if there are any overly long edges that we should subdivide
        for edge in list(face.edges):
            length = (edge.verts[0].co - edge.verts[1].co).magnitude
            cuts = round(length / targetlen) - 1
            if cuts >= 1:
                bmesh.ops.subdivide_edges(bm, edges = [edge], cuts = cuts)
        
        # Cast the face into 2D coordinate space
        poly = Polygon([v.co for v in face.verts])
        
        # Build a mapping from resulting 2D points back to the original vertices.
        # We need this when building 3D faces later.
        vertidx = dict((id(poly.points[i]), face.verts[i]) for i in range(len(face.verts)))
        
        print("Loading: %6.3f s" % (time.time() - start))
        start = time.time()
        
        # Generate vertices at each triangle grid intersection inside the polygon
        newpoints = []
        locationidx = {}
        xstepsize = targetlen
        ystepsize = math.sqrt(3)/2 * xstepsize
        ysteps = math.ceil((poly.max_y - poly.min_y) / ystepsize)
        xsteps = math.ceil((poly.max_x - poly.min_x) / xstepsize)
        offsety = ((poly.max_y - poly.min_y) - (ysteps * ystepsize)) / 2
        offsetx = ((poly.max_x - poly.min_x) - (xsteps * xstepsize)) / 2
        for yidx in range(ysteps):
            y = poly.min_y + offsety + (yidx + 0.5) * ystepsize
            for xidx in range(xsteps):
                mxidx = xidx * 2 + (1 - yidx % 2)
                x = poly.min_x + offsetx + mxidx * xstepsize / 2
                
                p = Vector([x, y])
                if poly.contains(p) and poly.edge_distance(p)[1] > xstepsize / 4:
                    locationidx[(mxidx, yidx)] = p
                    newpoints.append(p)
        
        # Index from point id() to point instance.
        # To make it possible to ignore duplicate edges, the items must be hashable, so
        # the code above uses id() of points.
        # Perhaps freezing the vectors could work also.
        allpoints = newpoints + poly.points
        pointidx = dict((id(p), p) for p in allpoints)
        orig_edges = {sortedge((id(p1), id(p2))) for p1, p2 in iter_tuples(poly.points, 2)}
        
        # Generate faces and edges for the regularly spaced inner area.
        regular_edges = set()
        regular_faces = set()
        for xidx, yidx in locationidx.keys():
            # Check the triangles starting in +Y and -Y directions from here
            # p1--p2
            #  \ /
            #   c
            #  / \
            # p3--p4
            
            p1 = locationidx.get((xidx-1,yidx-1))
            p2 = locationidx.get((xidx+1,yidx-1))
            c = locationidx.get((xidx,yidx))
            p3 = locationidx.get((xidx-1,yidx+1))
            p4 = locationidx.get((xidx+1,yidx+1))
            
            if p1 and p2 and c:
                regular_edges |= {sortedge((id(p1), id(p2))),
                                  sortedge((id(p2), id(c))),
                                  sortedge((id(c), id(p1)))}
                regular_faces |= {sortface(pointidx, (id(p1), id(p2), id(c)))}
            
            if p3 and p4 and c:
                regular_edges |= {sortedge((id(p3), id(p4))),
                                  sortedge((id(p4), id(c))),
                                  sortedge((id(c), id(p3)))}
                regular_faces |= {sortface(pointidx, (id(p3), id(p4), id(c)))}
        
        # Figure out which vertices lie outside the regular inner area
        edgecounts = dict((id(p),0) for p in newpoints)
        for p1, p2 in regular_edges:
            edgecounts[p1] += 1
            edgecounts[p2] += 1
        borderpoints = {p for p, c in edgecounts.items() if c < 6}
        regular_edges_on_border = {e for e in regular_edges
                                     if e[0] in borderpoints and e[1] in borderpoints}
        
        print("Regulars: %6.3f s" % (time.time() - start))
        start = time.time()
        
        # Connecting the border area vertices:
        # Generate edges between vertices that are close enough to each other
        border_edge_candidates = set()
        for p1 in [id(p) for p in poly.points]:
            for p2 in borderpoints:
                newedge = sortedge((p1,p2))
                if newedge not in regular_edges and newedge not in orig_edges:
                    if (pointidx[p1] - pointidx[p2]).magnitude <= targetlen * 2:
                        border_edge_candidates.add(newedge)
        
        # When edges intersect each other, keep only the shortest edge.
        discardededges = set()
        border_edge_candidates = list(border_edge_candidates)
        for i, e1 in enumerate(border_edge_candidates):
            for e2 in regular_edges_on_border:
                if e1[0] not in e2 and e1[1] not in e2:
                    points = [pointidx[e1[0]], pointidx[e1[1]], pointidx[e2[0]], pointidx[e2[1]]]
                    if geometry.intersect_line_line_2d(*points):
                        # Always favor the regular edges
                        discardededges.add(e1)
                        break
            
            if e1 not in discardededges:
                for e2 in border_edge_candidates[i+1:]:
                    if (e1[0] not in e2) and (e1[1] not in e2) and (e2 not in discardededges):
                        points = [pointidx[e1[0]], pointidx[e1[1]], pointidx[e2[0]], pointidx[e2[1]]]
                        if geometry.intersect_line_line_2d(*points):
                            # Keep shorter one
                            d1 = (pointidx[e1[0]] - pointidx[e1[1]]).magnitude
                            d2 = (pointidx[e2[0]] - pointidx[e2[1]]).magnitude
                            if d1 < d2:
                                discardededges.add(e2)
                            else:
                                discardededges.add(e1)
        border_edges = {e for e in border_edge_candidates if e not in discardededges}

        print("Border edges: %6.3f s" % (time.time() - start))
        start = time.time()

        # Build a lookup from point id() to edge tuples
        # It allows us to find all edges that originate from a given point.
        alledges = regular_edges | border_edges | orig_edges
        edgelookup = dict((id(p),set()) for p in allpoints)
        for e in alledges:
            for p in e:
                edgelookup[p].add(e)

        # Take one edge at a time, and then try to find two other edges that
        # share points. Once such a triplet is found, make a triangular face
        # out of it.
        border_faces = set()
        for e1 in border_edges:
            e2_candidates = [e for e in edgelookup[e1[0]] if e1 != e]
            e3_candidates = [e for e in edgelookup[e1[1]] if e1 != e]
            for e2 in e2_candidates:
                p3 = e2[0] if e2[0] != e1[0] else e2[1]
                e3 = [e for e in e3_candidates if p3 in e]
                if e3:
                    border_faces.add(sortface(pointidx, (e1[0], e1[1], p3)))
    
        print("Border faces: %6.3f s" % (time.time() - start))
        start = time.time()
    
        # Add all new points to the bmesh as 3D vertices.
        # Update the vertidx mapping also.
        for p in newpoints:
            vertidx[id(p)] = bm.verts.new(poly.plane.to_3d(p))
        
        if regular_faces or border_faces:
            # Remove the original selected face
            bmesh.ops.delete(bm, geom=[bm.faces.active], context=3)
        
        if self.only_edges:
            # For debugging: show the edges instead of faces
            for p1, p2 in border_edges | regular_edges:
                print("Adding " + str((p1,p2)))
                bm.edges.new((vertidx[p1], vertidx[p2]))
        else:
            # Add all the new faces by looking up 3d vertices using the map.
            for p1, p2, p3 in (regular_faces | border_faces):
                f = bm.faces.new((vertidx[p1], vertidx[p2], vertidx[p3]))
                f.normal_update()
        
        bmesh.update_edit_mesh(obj.data)
        
        print("Storage: %6.3f s" % (time.time() - start))
        start = time.time()
        return {'FINISHED'}

def register():
    bpy.utils.register_class(TriangleFill)

def unregister():
    bpy.utils.unregister_class(TriangleFill)

if __name__ == "__main__":
    register()


