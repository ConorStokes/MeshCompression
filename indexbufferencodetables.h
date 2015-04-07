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

// NOTE - these tables are in this file for readability / navagability and should *only* be included in
// indexbuffercompression.cpp for index buffer compression and meshcompression.cpp for mesh compression


#pragma once

// Individual vertex type classifications.
enum VertexClassification
{
    NEW_VERTEX = 0,
    CACHED_VERTEX = 1,
    FREE_VERTEX = 2
};

// Individual case for handling a combination of vertice classifications.
struct VertexCompressionCase
{
    IndexBufferTriangleCodes code;
    uint32_t vertexOrder[ 3 ];
};

// This is a table for looking up the appropriate code and rotation for a set of vertex classifications.
static const VertexCompressionCase CompressionCase[3][3][3] =
{
    { // new 
        { // new new
            { // new new new 
                IB_NEW_NEW_NEW, { 0, 1, 2 }
            },
            { // new new cached
                IB_NEW_NEW_CACHED, { 0, 1, 2 }
            },
            { // new new free
                IB_NEW_NEW_FREE, { 0, 1, 2 }
            }
        },
        { // new cached
            { // new cached new
                IB_NEW_NEW_CACHED, { 2, 0, 1 }
            },
            {  // new cached cached
                IB_NEW_CACHED_CACHED, { 0, 1, 2 }
            },
            { // new cached free
                IB_NEW_CACHED_FREE, { 0, 1, 2 }
            }
        },
        { // new free
            { // new free new
                IB_NEW_NEW_FREE, { 2, 0, 1 }
            },
            { // new free cached
                IB_NEW_FREE_CACHED, { 0, 1, 2 }
            },
            { // new free free
                IB_NEW_FREE_FREE, { 0, 1, 2 }
            }
        }
    },
    { // cached
        { // cached new 
            { // cached new new
                IB_NEW_NEW_CACHED, { 1, 2, 0 }
            },
            { // cached new cached
                IB_NEW_CACHED_CACHED, { 1, 2, 0 }
            },
            { // cached new free
                IB_NEW_FREE_CACHED, { 1, 2, 0 }
            }
        },
        { // cached cached
            { // cached cached new
                IB_NEW_CACHED_CACHED, { 2, 0, 1 }
            },
            { // cached cached cached
                IB_CACHED_CACHED_CACHED, { 0, 1, 2 }
            },
            { // cached cached free
                IB_CACHED_CACHED_FREE, { 0, 1, 2 }
            }
        },
        { // cached free
            { // cached free new
                IB_NEW_CACHED_FREE, { 2, 0, 1 }
            },
            { // cached free cached
                IB_CACHED_CACHED_FREE, { 2, 0, 1 }
            },
            { // cached free free 
                IB_CACHED_FREE_FREE, { 0, 1, 2 }
            }
        }
    },
    { // free
        { // free new
            { // free new new
                IB_NEW_NEW_FREE, { 1, 2, 0 }
            },
            { // free new cached
                IB_NEW_CACHED_FREE, { 1, 2, 0 }
            },
            { // free new free
                IB_NEW_FREE_FREE, { 1, 2, 0 }
            }
        },
        { // free cached
            { // free cached new
                IB_NEW_FREE_CACHED, { 2, 0, 1 }
            },
            { // free cached cached
                IB_CACHED_CACHED_FREE, { 1, 2, 0 }
            },
            { // free cached free
                IB_CACHED_FREE_FREE, { 1, 2, 0 }
            }
        },
        { // free free
            { // free free new
                IB_NEW_FREE_FREE, { 2, 0, 1 }
            },
            { // free free cached
                IB_CACHED_FREE_FREE, { 2, 0, 1 }
            },
            { // free free free
                IB_FREE_FREE_FREE, { 0, 1, 2 }
            }
        }
    }
};

// Prefix code table used for encoding edge bits
static const PrefixCode EdgePrefixCodes[] =
{
    { 1, 2 },
    { 2, 2 },
    { 0, 3 },
    { 15, 4 },
    { 11, 4 },
    { 3, 4 },
    { 7, 5 },
    { 28, 5 },
    { 20, 5 },
    { 55, 6 },
    { 12, 6 },
    { 36, 6 },
    { 23, 7 },
    { 44, 7 },
    { 215, 8 },
    { 87, 8 },
    { 196, 8 },
    { 132, 8 },
    { 236, 9 },
    { 364, 9 },
    { 324, 9 },
    { 68, 9 },
    { 1004, 10 },
    { 492, 10 },
    { 108, 10 },
    { 772, 10 },
    { 516, 10 },
    { 4, 10 },
    { 1644, 11 },
    { 620, 11 },
    { 1284, 11 },
    { 260, 11 }
};

// Prefix code table used for vertices
static const PrefixCode CachedVertexPrefixCodes[] =
{
	{ 215, 8 },
	{ 0, 1 },
	{ 5, 3 },
	{ 3, 4 },
	{ 15, 5 },
	{ 11, 5 },
	{ 9, 5 },
	{ 1, 5 },
	{ 55, 6 },
	{ 39, 6 },
	{ 27, 6 },
	{ 25, 6 },
	{ 17, 6 },
	{ 63, 7 },
	{ 31, 7 },
	{ 23, 7 },
	{ 7, 7 },
	{ 59, 7 },
	{ 121, 7 },
	{ 113, 7 },
	{ 49, 7 },
	{ 255, 8 },
	{ 127, 8 },
	{ 223, 8 },
	{ 95, 8 },
	{ 87, 8 },
	{ 199, 8 },
	{ 71, 8 },
	{ 251, 8 },
	{ 123, 8 },
	{ 185, 8 },
	{ 57, 8 }
};

// Prefix code table used for triangles
static const PrefixCode TrianglePrefixCodes[] =
{
	{ 0, 1 },
	{ 3, 2 },
	{ 5, 3 },
	{ 49, 7 },
	{ 33, 7 },
	{ 81, 7 },
	{ 9, 5 },
	{ 113, 7 },
	{ 57, 7 },
	{ 25, 6 },
	{ 121, 7 },
	{ 17, 7 },
	{ 1, 6 },
	{ 97, 7 }
};

// Constant value for vertices that don't get mapped in the vertex re-map.
static const uint32_t VERTEX_NOT_MAPPED = 0xFFFFFFFF;