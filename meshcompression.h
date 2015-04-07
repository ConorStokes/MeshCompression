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
#ifndef MESH_COMPRESSION_H__
#define MESH_COMPRESSION_H__
#pragma once

#include <stdint.h>

class WriteBitstream;

// Compress an index buffer, writing the results out to a bitstream and providing a vertex remapping (which will be in pre-transform cache optimised
// order).
//
// Recommended maximum range for vertex attributes it -2^29 to 2^29 - 1 
//
// Parameters: 
//     [in]  triangles           - A typical triangle list index buffer (3 indices to vertices per triangle). 16 bit indices.
//     [in]  triangle count      - The number of triangles to process.
//     [out] vertexRemap         - This will be populated with re-mappings that map old vertices to new vertex locations (a new ordering),
//                                 where indexing with the old vertex index will get you the new one. Vertices that are unused will 
//                                 be mapped to 0xFFFFFFFF.
//                                 You should re-order the vertices and removed unused ones based on the vertex remap, instead of storing
//                                 the remap. 
//                                 It should be allocated as a with at least vertexCount entries.
//     [in] vertexCount          - The number of vertices in the mesh. This should be less than 0xFFFFFFFF/2^32 - 1.
//     [in] vertexAttributeCount - The number of attributes for each vertice in the mesh. 
//     [in] vertexAttributes     - The vertex attributes (the attributes for each vertex are packed together, so there are vertexCount * vertexAttributeCount entries.
//     [in] output               - The stream that the compressed data will be written to. Note that we will not flush/finish the stream
//                                 in case something else is going to be written after, so WriteBitstream::Finish will need to be called after this.
void CompressMesh( 
	const uint16_t* triangles, 
	uint32_t triangleCount, 
	uint32_t* vertexRemap, 
	uint32_t vertexCount, 
	uint32_t vertexAttributeCount, 
	const int32_t* vertexAttributes, 
	WriteBitstream& output );

// Same as above but 32bit indices.
void CompressMesh(
	const uint32_t* triangles,
	uint32_t triangleCount,
	uint32_t* vertexRemap,
	uint32_t vertexCount,
	uint32_t vertexAttributeCount,
	const int32_t* vertexAttributes,
	WriteBitstream& output );

// Same as above but 16bit indices and 16 bit attributes.
// Recommended maximum range for vertex attributes it -2^14 to 2^14 - 1, if you wish to use the 16bit decoder.
void CompressMesh(
	const uint16_t* triangles,
	uint32_t triangleCount,
	uint32_t* vertexRemap,
	uint32_t vertexCount,
	uint32_t vertexAttributeCount,
	const int16_t* vertexAttributes,
	WriteBitstream& output );

// Same as above but 32bit indices and 16 bit attributes.
// Recommended maximum range for vertex attributes it -2^14 to 2^14 - 1, if you wish to use the 16bit decoder.
void CompressMesh(
	const uint32_t* triangles,
	uint32_t triangleCount,
	uint32_t* vertexRemap,
	uint32_t vertexCount,
	uint32_t vertexAttributeCount,
	const int16_t* vertexAttributes,
	WriteBitstream& output );

#endif