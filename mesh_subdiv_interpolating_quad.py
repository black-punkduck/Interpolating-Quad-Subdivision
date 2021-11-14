
bl_info = {
    "name": "Interpolating Subdivision for Quads",
    "author": "black-punkduck",
    "version": (0, 9, 0),
    "blender": (2, 82, 0),
    "location": "Object > Transform > Interpolating Subdivision for Quads",
    "description": "Applies Quad Subdivision with Kobbelt algorithm",
    "category": "Object",
}

import bpy
import bmesh
import mathutils
from bpy.types import Operator
from bpy.utils import register_class, unregister_class

class KobbeltSubdiv:
    def __init__(self, obj):
        self.obj = obj
        self.data = obj.data
        self.bm = bmesh.new()

        # convert the current mesh to a bmesh (must be in edit mode)
        bpy.ops.object.mode_set(mode='EDIT')
        self.bm.from_mesh(self.data)
        bpy.ops.object.mode_set(mode='OBJECT')  # return to object mode

        self.vertEdges = {}        # contains all edges connected to a vertex
        self.vertPolygons = {}     # contains all polygons connected to a vertex

        self.edgePolygons = {}     # contains all polygons connected to an edge

        self.polygonEdges = {}     # contains all edges a polygon is using
        self.polygonNeighbors = {} # contains all neighbor-polygons/faces
        self.polygonUVVerts = {}   # contains the UV-Vertices of all polygons for later use

        self.newEdgeMCoord = {}  # contains all new mid-vertex coordinates of the subdivided edges
        self.newEdgeMIndex = []  # contains all indices of the new vertices on the edges
        self.newPolyMCoord = {}  # contains all new mis-poly coordinatese of the subdivided quad
        self.newPolyMIndex = []  # contains all indices of the new vertices on the edges

        self.newEdges = []  # creates all elements for drawing new Edges

        self.uvlayer  = self.data.uv_layers.active


    def closemesh(self):
        self.bm.to_mesh(self.data)
        self.bm.free()

    def calculatehelper(self):
        mesh = self.data

        for vertex in mesh.vertices:
            self.vertEdges[vertex.index] = []
            self.vertPolygons[vertex.index] = []

        # supply an index of all polygons connected to a vertex
        # and check if number of vertices is always 4

        for polygon in mesh.polygons:
            if len(polygon.vertices) != 4:
                return (False)
            for vertex in polygon.vertices:
                self.vertPolygons[vertex].append(polygon)

        # supply an index of all edges connected to a vertex

        for edge in mesh.edges:
            for vertex in edge.vertices:
                self.vertEdges[vertex].append(edge)


        # supply an index of all edges connected to a polygon
        #
        for edge in mesh.edges:
            self.edgePolygons[edge.index] = []

        self.newEdgeMIndex = [None] * len(mesh.edges)
        for polygon in mesh.polygons:
            self.polygonEdges[polygon.index] = []

        for polygon in mesh.polygons:
            for vertex in polygon.vertices:
                for edge in self.vertEdges[vertex]:
                    v0 = edge.vertices[0]
                    v1 = edge.vertices[1]
                    if (v0 in polygon.vertices) and (v1 in polygon.vertices):
                        if polygon not in self.edgePolygons[edge.index]:
                            self.edgePolygons[edge.index].append(polygon)
                        if edge not in self.polygonEdges[polygon.index]:
                            self.polygonEdges[polygon.index].append(edge)

        # evaluate polygon neighbors (or faces neighbors)
        #
        # an entry containing edge number and polygon number is created

        self.newPolyMIndex = [None] * len(mesh.polygons)

        for polygon in mesh.polygons:
            self.polygonNeighbors[polygon.index] = {}
            self.polygonUVVerts[polygon.index] = []

        for polygon in mesh.polygons:
            for edge in self.polygonEdges[polygon.index]:
                for neighbor in self.edgePolygons[edge.index]:
                    if neighbor != polygon:                 # do not add yourself :-)
                        self.polygonNeighbors[polygon.index][edge.index] = neighbor

        #
        # keep the old polygon UV vertices
        #
        if self.uvlayer is not None:
            for polygon in mesh.polygons:
                for vert_idx, loop_idx in zip(polygon.vertices, polygon.loop_indices):
                    self.polygonUVVerts[polygon.index].append(self.uvlayer.data[loop_idx].uv)

        return (True)

    def getNeighborVertex(self, p_id, edge_id):
        """
        method to find a neighbor, which is not part of the polygons linked
        to the same edge
        """
        for edge in self.vertEdges[p_id]:
            #print("Testing Edge: " + str(edge.index))
            found = 0
            for poly in self.edgePolygons[edge_id]:
                if poly in self.edgePolygons[edge.index]:
                    found = 1
                    continue
            if found == 1:
                continue
            # v0 oder v1
            res_id = edge.vertices[0] if (p_id != edge.vertices[0]) else edge.vertices[1]
            return (res_id)

    def getNeighborSubdiv(self, polyneighbors, pid):
        """
        method to find edge on the opposite side of a neighbor polygon
        """
        prev = polyneighbors[pid].index
        #print ("Neighbor: " + str(prev))
        pe = [edge.index for edge in self.polygonEdges[prev]]
        #print (pe)
        idx = (pe.index(pid) + 2) & 3
        #print ("Index is" + str(idx))
        p = self.newEdgeMCoord[pe[idx]]
        return (p)

    def createMidEdgeVerts(self):
        """
        calculate new vertices on old edges
        Append the verts to an array to draw later
        """
        print ("Create mid vertices for edges")
        mesh = self.data
        for edge in mesh.edges:
            if self.newEdgeMCoord[edge.index] is not None:
                num = self.bm.verts.new(self.newEdgeMCoord[edge.index]) 
                self.bm.verts.index_update()
                self.newEdges.append((edge.vertices[0], num.index))
                self.newEdges.append((num.index, edge.vertices[1]))
                self.newEdgeMIndex[edge.index] = num.index

    def createMidPolygonVerts(self):
        """
        calculate new vertices to connect to mid polygon
        Append the verts to an array to draw later
        """
        print ("Create mid vertices for quads")
        mesh = self.data
        for poly in mesh.polygons:
            if self.newPolyMCoord[poly.index] is not None:
                num = self.bm.verts.new(self.newPolyMCoord[poly.index]) 
                self.bm.verts.index_update()
                self.newPolyMIndex[poly.index] = num.index
                for edge in self.polygonEdges[poly.index]:
                    self.newEdges.append((self.newEdgeMIndex[edge.index], num.index))



    def splitedge(self):
        mesh = self.data
        bm   = self.bm
        #oldverts = bm.verts
        oldverts = mesh.vertices

        for edge in mesh.edges:
            print ("Edge: " + str(edge.index))
            # print (self.edgePolygons[edge.index])   # anliegende Polygone

            p2_id = edge.vertices[0]
            p2_polenum = len(self.vertEdges[p2_id])
            p2 = oldverts[p2_id].co

            p3_id = edge.vertices[1]
            p3_polenum = len(self.vertEdges[p3_id])
            p3 = oldverts[p3_id].co

            valid = 1
            p1_id = -1
            p1 = mathutils.Vector((0.0,0.0,0.0))
            p4_id = -1
            p4 = mathutils.Vector((0.0,0.0,0.0))

            if p2_polenum >= 4:
                p1_id = self.getNeighborVertex(p2_id, edge.index)
                p1 = oldverts[p1_id].co

            elif p2_polenum == 2:
                p1 = (2.0 * p2) - p3
            else:
                print ("Polenum = " + str(p2_polenum))
                p1 = (2.0 * p2) - p3
                #valid = 0

            if p3_polenum >= 4:
                p4_id = self.getNeighborVertex(p3_id, edge.index)
                p4 = oldverts[p4_id].co

            elif p3_polenum == 2:
                p4 = (2.0 * p3) - p2
            else:
                print ("Polenum = " + str(p3_polenum))
                p4 = (2.0 * p3) - p2
                #valid = 0


            print("p1_id = " + str(p1_id) + " " + str(p1))
            print("p2_id = " + str(p2_id) + " " + str(p2))
            print("p3_id = " + str(p3_id) + " " + str(p3))
            print("p4_id = " + str(p4_id) + " " + str(p4))
            if valid == 0:
                self.newEdgeMCoord[edge.index] = None
            else:
                newvert = 0.5625 * (p2 + p3) - 0.0625 * (p1 + p4)
                self.newEdgeMCoord[edge.index] = newvert
        
        for poly in mesh.polygons:
            print ("Polygon: " + str(poly.index))
            pn = self.polygonNeighbors[poly.index]
            pe = [edge.index for edge in self.polygonEdges[poly.index]]
            # neigbor for index 0 and 2, same for 1 and 3
            p1 = mathutils.Vector((0.0,0.0,0.0))
            p2 = mathutils.Vector((0.0,0.0,0.0))
            p3 = mathutils.Vector((0.0,0.0,0.0))
            p4 = mathutils.Vector((0.0,0.0,0.0))
            q1 = mathutils.Vector((0.0,0.0,0.0))
            q2 = mathutils.Vector((0.0,0.0,0.0))
            q3 = mathutils.Vector((0.0,0.0,0.0))
            q4 = mathutils.Vector((0.0,0.0,0.0))
            valid = 1
            if pe[0] in pn:
                p2 = self.newEdgeMCoord[pe[0]]
                p1 = self.getNeighborSubdiv(pn, pe[0])
                print ("P2 = " + str(p2) + ", P1 = " + str(p1))
            else:
                print ("No neighbor of e0")
                valid = 0

            if pe[1] in pn:
                q2 = self.newEdgeMCoord[pe[1]]
                q1 = self.getNeighborSubdiv(pn, pe[1])
                print ("Q2 = " + str(q2) + ", Q1 = " + str(q1))
            else:
                print ("No neighbor of e1")
                valid = 0

            if pe[2] in pn:
                p3 = self.newEdgeMCoord[pe[2]]
                p4 = self.getNeighborSubdiv(pn, pe[2])
                print ("P3 = " + str(p3) + ", P4 = " + str(p4))
            else:
                print ("No neighbor of e2")
                valid = 0

            if pe[3] in pn:
                q3 = self.newEdgeMCoord[pe[3]]
                q4 = self.getNeighborSubdiv(pn, pe[3])
                print ("Q3 = " + str(q3) + ", Q4 = " + str(q4))
            else:
                print ("No neighbor of e3")
                valid = 0


            if valid == 1:
                p = 0.5625 * (p2 + p3) - 0.0625 * (p1 + p4)
                q = 0.5625 * (q2 + q3) - 0.0625 * (q1 + q4)
                self.newPolyMCoord[poly.index] = 0.5 * (p + q)
            else:
                self.newPolyMCoord[poly.index] = None

        self.createMidEdgeVerts()
        self.createMidPolygonVerts()

        print ("Delete old edges")
        #bmesh.ops.delete(bm, geom=[edge for edge in bm.edges], context=4)
        bmesh.ops.delete(bm, geom=[edge for edge in bm.edges], context='EDGES_FACES')

        print ("Create new edges")
        bm.verts.ensure_lookup_table()
        for v in self.newEdges:
            bm.edges.new((bm.verts[v[0]], bm.verts[v[1]]))
         
        bm.verts.ensure_lookup_table()

        print ("Create new faces")

        uvl = bm.loops.layers.uv[0]

        for poly in mesh.polygons:
            # print (poly.index)
            midpoly = self.newPolyMIndex[poly.index]
            old_uvs = self.polygonUVVerts[poly.index]
            if midpoly is not None:
                #
                # each edge draws one polygon
                #
                miduvp = (old_uvs[0] + old_uvs[1] + old_uvs[2] + old_uvs[3]) / 4
                for i,v1 in enumerate(poly.vertices):
                    i2 = (i+1) & 3
                    i3 = (i+3) & 3
                    v2 = poly.vertices[i2]
                    v3 = poly.vertices[i3]
                    midedge = None
                    midedge2 = None
                    miduv1 = mathutils.Vector((0.0,0.0))
                    miduv2 = mathutils.Vector((0.0,0.0))
                    for edge in self.polygonEdges[poly.index]:
                        if v1 in edge.vertices and v2 in edge.vertices:
                            midedge = self.newEdgeMIndex[edge.index]
                            miduv1 = (old_uvs[i] + old_uvs[i2]) / 2
                            break
                    for edge in self.polygonEdges[poly.index]:
                        if v1 in edge.vertices and v3 in edge.vertices:
                            midedge2 = self.newEdgeMIndex[edge.index]
                            miduv2 = (old_uvs[i] + old_uvs[i3]) / 2
                            break

                    try:
                        fid = bm.faces.new([bm.verts[v1], bm.verts[midedge], bm.verts[midpoly], bm.verts[midedge2]])
                    except:
                        print ("Ignore that guy")

                    fid.loops[0][uvl].uv = old_uvs[i]
                    fid.loops[1][uvl].uv = miduv1
                    fid.loops[2][uvl].uv = miduvp
                    fid.loops[3][uvl].uv = miduv2
    


class SubDivKobbelt(Operator):
    bl_idname = "subdiv.kobbelt"
    bl_label = "Kobbelt algorithm"
    bl_options = {'REGISTER', 'UNDO'}


    @classmethod
    def poll(cls, context):
        obj = context.active_object
        return (obj and obj.type == 'MESH')

    def execute(self, context):
        
        print ("Kobbelt called\n")

        subdiv = KobbeltSubdiv(context.active_object)
        if subdiv.calculatehelper() is False:
            self.report({'ERROR'}, "Mesh is not quad.")
            return {'FINISHED'}
        subdiv.splitedge()
        subdiv.closemesh()
        return {'FINISHED'}

def menu_func(self, context):
    self.layout.operator("subdiv.kobbelt", text="Kobbelt Subdivsion")

def register():
    bpy.utils.register_class(SubDivKobbelt)
    bpy.types.VIEW3D_MT_object.append(menu_func)


def unregister():
    bpy.utils.unregister_class(SubDivKobbelt)
    bpy.types.VIEW3D_MT_object.remove(menu_func)


if __name__ == "__main__":
    register()

