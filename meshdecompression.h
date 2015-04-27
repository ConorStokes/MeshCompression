/*
Copyright (c) 2015, Conor Stokes
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

1. Redistributions of source code must retain the above copyright notice, this
list of conditions and the following disclaimer.
2. Redistributions in binary form must reproduce the above copyright notice,
this list of conditions and the following disclaimer in the documentation
and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND
ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED
WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR
ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
(INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef MESH_DECOMPRESSION_H__
#define MESH_DECOMPRESSION_H__
#pragma once

#include <stdint.h>
#include "readbitstream.h"

// Decompress a triangle mesh, consisting of a set of vertices, referenced by a list of triangles (indices)
// The vertices consist of a set of vertex attributes, all are int32_ts.
// Note, vertex attributes come out in their quantitized form, as they went in to compress mesh
// All vertex attributes use delta coding using either a parallelogram predictor (for edge cache hits)
// or another vertex in the triangle (except for NEW NEW NEW cases, where the first vertex is encoded in absolute terms).
// Recommended maximum range for attributes is -2^29 to 2^29 - 1. 
// Parameters: 
//     [out] triangles            - Triangle list index buffer (3 indices to vertices per triangle), output from the decompression - 16bit indices
//     [in]  triangleCount        - The number of triangles to decompress.
//     [in]  vertexAttributeCount - The number of attributes per vertex
//     [out] vertexAttributes     - The decompressed vertex attributes.
//     [in]  input                - The bit stream that the compressed data will be read from.
void DecompressMesh(
    uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int32_t* vertexAttributes,
    ReadBitstream& input );

// Same as above but 16 bit indices.
void DecompressMesh(
    uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int32_t* vertexAttributes,
    ReadBitstream& input );

// Same as above but 32 bit indices and 16 bit vertex attributes. 
void DecompressMesh(
    uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int16_t* vertexAttributes,
    ReadBitstream& input );

// Same as above but 32 bit indices.
void DecompressMesh(
    uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int16_t* vertexAttributes,
    ReadBitstream& input );


#endif // -- MESH_DECOMPRESSION_H__
