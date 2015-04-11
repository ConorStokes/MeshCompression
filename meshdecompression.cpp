/*
Copyright (c) 2014-2015, Conor Stokes
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
#include "indexbufferdecompression.h"
#include "readbitstream.h"
#include "indexcompressionconstants.h"
#include "indexbuffercompressionformat.h"
#include "meshcompressionconstants.h"
#include <assert.h>

static const uint32_t EDGE_MAX_CODE_LENGTH     = 11;
static const uint32_t VERTEX_MAX_CODE_LENGTH   = 8;
static const uint32_t TRIANGLE_MAX_CODE_LENGTH = 7;

#ifdef _MSC_VER
#define MDC_INLINE __forceinline
#else
#define MDC_INLINE inline
#endif 

#include "indexbufferdecodetables.h"

template < typename AttributeType >
AttributeType DecodeExponentialGolomb( ReadBitstream& input, const uint32_t k, /* out */ uint32_t& bitSize )
{
    return 0;
}

template <>
MDC_INLINE int16_t DecodeExponentialGolomb< int16_t >( ReadBitstream& input, const uint32_t k, /* out */ uint32_t& bitSize )
{
    return input.DecodeExponentialGolomb16( k, bitSize );
}

template <>
MDC_INLINE int32_t DecodeExponentialGolomb< int32_t >( ReadBitstream& input, const uint32_t k, /* out */ uint32_t& bitSize )
{
    return input.DecodeExponentialGolomb( k, bitSize );
}

// Decompress triangle codes using prefix coding based on static tables.
template <typename IndiceType, typename AttributeType>
void DecompressMeshPrefix( 
    IndiceType* triangles,
    uint32_t triangleCount, 
    uint32_t vertexAttributeCount,
    AttributeType* vertexAttributes,
    ReadBitstream& input )
{
    EdgeTriangle edgeFifo[ EDGE_FIFO_SIZE ];
    uint32_t     vertexFifo[ VERTEX_FIFO_SIZE ];

    uint32_t          edgesRead    = 0;
    uint32_t          verticesRead = 0;
    uint32_t          newVertices  = 0;
    const IndiceType* triangleEnd  = triangles + ( triangleCount * 3 );
    AttributeType*    newVertex    = vertexAttributes;

    // array of exponential moving average values to estimate optimal k for exp golomb codes
    // note we use 16/16 unsigned fixed point.
    uint32_t kArray[ 64 ];

    for ( uint32_t where = 0; where < vertexAttributeCount; ++where )
    {
        // prime the array of ks for exp golomb with an average bitsize of 4
        // note that k is 
        kArray[ where ] = 4 << 16;
    }

    // iterate through the triangles
    for ( IndiceType* triangle = triangles; triangle < triangleEnd; triangle += 3 )
    {
        IndexBufferTriangleCodes code = static_cast< IndexBufferTriangleCodes >( input.Decode( TriangleDecoding, TRIANGLE_MAX_CODE_LENGTH ) );

        switch ( code )
        {
        case IB_EDGE_NEW:
        {
            uint32_t            edgeFifoIndex = input.Decode( EdgeDecoding, EDGE_MAX_CODE_LENGTH );
            const EdgeTriangle& edge          = edgeFifo[ ( ( edgesRead - 1 ) - edgeFifoIndex ) & EDGE_FIFO_MASK ];

            triangle[ 0 ]                               = static_cast< IndiceType >( edge.second );
            triangle[ 1 ]                               = static_cast< IndiceType >( edge.first );

            vertexFifo[ verticesRead & EDGE_FIFO_MASK ] =
            triangle[ 2 ]                               = static_cast< IndiceType >( newVertices );

            const AttributeType* adjacent1Attribute = vertexAttributes + ( edge.first * vertexAttributeCount );
            const AttributeType* adjacent2Attribute = vertexAttributes + ( edge.second * vertexAttributeCount );
            const AttributeType* opposingAttribute  = vertexAttributes + ( edge.third * vertexAttributeCount );
            AttributeType*       endAttributes      = newVertex + vertexAttributeCount;
            uint32_t*            k                  = kArray;

            for ( ; newVertex < endAttributes; ++adjacent1Attribute, ++adjacent2Attribute, ++opposingAttribute, ++newVertex, ++k )
            {
                AttributeType predicted = *adjacent2Attribute + ( *adjacent1Attribute  - *opposingAttribute );
                uint32_t kEstimate;

                *newVertex = predicted + DecodeExponentialGolomb< AttributeType >( input, *k >> 16, kEstimate );

                // fixed point exponential moving average with alpha 0.0625 (equivalent to N being 31)
                *k = ( *k * 15 + ( kEstimate << 16 ) ) >> 4;
            }

            ++newVertices;
            ++verticesRead;

            break;
        }

        case IB_EDGE_CACHED:
        {
            uint32_t            edgeFifoIndex   = input.Decode( EdgeDecoding, EDGE_MAX_CODE_LENGTH );
            uint32_t            vertexFifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            const EdgeTriangle& edge            = edgeFifo[ ( ( edgesRead - 1 ) - edgeFifoIndex ) & EDGE_FIFO_MASK ];

            triangle[ 0 ] = static_cast< IndiceType >( edge.second );
            triangle[ 1 ] = static_cast< IndiceType >( edge.first );
            triangle[ 2 ] = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertexFifoIndex ) & VERTEX_FIFO_MASK ] );

            break;
        }
        case IB_EDGE_FREE:
        {
            uint32_t            edgeFifoIndex   = input.Decode( EdgeDecoding, EDGE_MAX_CODE_LENGTH );
            uint32_t            relativeVertex  = input.ReadVInt();
            const EdgeTriangle& edge            = edgeFifo[ ( ( edgesRead - 1 ) - edgeFifoIndex ) & EDGE_FIFO_MASK ];

            triangle[ 0 ]                                 = static_cast< IndiceType >( edge.second );
            triangle[ 1 ]                                 = static_cast< IndiceType >( edge.first );

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                 = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex );

            ++verticesRead;

            break;
        }
        case IB_NEW_NEW_NEW:
        {
            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( newVertices + 1 );
            vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( newVertices + 2 );
            
            const uint32_t* k        = kArray;

            AttributeType* vert0End = newVertex + vertexAttributeCount;

            for ( ; newVertex < vert0End; ++newVertex, ++k )
            {
                uint32_t dummy;
                AttributeType readVert0 = DecodeExponentialGolomb< AttributeType >( input, EXP_GOLOMB_FIRST_NEW_K, dummy );

                *newVertex                                    = readVert0;
                *( newVertex + vertexAttributeCount )         = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert0;
                *( newVertex + ( 2 * vertexAttributeCount ) ) = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert0;
            }
            
            newVertex += 2 * vertexAttributeCount;

            newVertices  += 3;
            verticesRead += 3;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_NEW_CACHED:
        {
            uint32_t vertexFifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );

            triangle[ 2 ]                                         = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertexFifoIndex ) & VERTEX_FIFO_MASK ] );
            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( newVertices + 1 );

            const AttributeType* vert2    = vertexAttributes + ( triangle[ 2 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert2, ++k )
            {
                AttributeType readVert2  = *vert2;
                uint32_t dummy;

                *newVertex                           = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert2;
                *(newVertex + vertexAttributeCount ) = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert2;
            }

            newVertex += vertexAttributeCount;

            verticesRead += 2;
            newVertices  += 2;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_NEW_FREE:
        {
            uint32_t relativeVertex = input.ReadVInt();

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( newVertices + 1 );
            vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex );

            const AttributeType* vert2    = vertexAttributes + ( triangle[ 2 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert2, ++k )
            {
                AttributeType readVert2 = *vert2;
                uint32_t dummy;

                *newVertex                           = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert2;
                *(newVertex + vertexAttributeCount ) = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + readVert2;
            }

            newVertex += vertexAttributeCount;

            newVertices  += 2;
            verticesRead += 3;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_CACHED_CACHED:
        {
            uint32_t vertex1FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t vertex2FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );

            triangle[ 1 ]                                 = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex1FifoIndex ) & VERTEX_FIFO_MASK ] );
            triangle[ 2 ]                                 = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex2FifoIndex ) & VERTEX_FIFO_MASK ] );
            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] =
            triangle[ 0 ]                                 = static_cast< IndiceType >( newVertices );

            const AttributeType* vert1    = vertexAttributes + ( triangle[ 1 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert1, ++k )
            {
                uint32_t dummy;

                *newVertex = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + *vert1;
            }
            
            ++verticesRead;
            ++newVertices;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_CACHED_FREE:
        {
            uint32_t vertexFifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t relativeVertex  = input.ReadVInt();

            triangle[ 1 ]                                         = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertexFifoIndex ) & VERTEX_FIFO_MASK ] );
            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex );		

            const AttributeType* vert1    = vertexAttributes + ( triangle[ 1 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert1, ++k )
            {
                uint32_t dummy;

                *newVertex = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + *vert1;
            }
            
            verticesRead += 2;
            ++newVertices;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_FREE_CACHED:
        {
            uint32_t relativeVertex  = input.ReadVInt();
            uint32_t vertexFifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );

            triangle[ 2 ]                                         = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertexFifoIndex ) & VERTEX_FIFO_MASK ] );
            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex );				

            const AttributeType* vert2    = vertexAttributes + ( triangle[ 2 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert2, ++k )
            {
                uint32_t dummy;

                *newVertex = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + *vert2;
            }
            
            verticesRead += 2;
            ++newVertices;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_NEW_FREE_FREE:
        {
            uint32_t relativeVertex1  = input.ReadVInt();
            uint32_t relativeVertex2  = input.ReadVInt();

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( newVertices );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex1 );
            vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex2 );
            
            const AttributeType* vert1    = vertexAttributes + ( triangle[ 1 ] * vertexAttributeCount );
            const AttributeType* vert0End = newVertex + vertexAttributeCount;
            const uint32_t*      k        = kArray;

            for ( ; newVertex < vert0End; ++newVertex, ++vert1, ++k )
            {
                uint32_t dummy;

                *newVertex = DecodeExponentialGolomb< AttributeType >( input, ( *k >> 16 ), dummy ) + *vert1;
            }
                        
            verticesRead += 3;
            ++newVertices;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;

            break;
        }
        case IB_CACHED_CACHED_CACHED:
        {
            uint32_t vertex0FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t vertex1FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t vertex2FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );

            triangle[ 0 ] = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex0FifoIndex ) & VERTEX_FIFO_MASK ] );
            triangle[ 1 ] = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex1FifoIndex ) & VERTEX_FIFO_MASK ] );
            triangle[ 2 ] = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex2FifoIndex ) & VERTEX_FIFO_MASK ] );

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;
            break;
        }
        case IB_CACHED_CACHED_FREE:
        {
            uint32_t vertex0FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t vertex1FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t relativeVertex2  = input.ReadVInt();

            triangle[ 0 ]                                 = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex0FifoIndex ) & VERTEX_FIFO_MASK ] );
            triangle[ 1 ]                                 = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex1FifoIndex ) & VERTEX_FIFO_MASK ] );

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                 = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex2 );

            ++verticesRead;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;

            break;
        }
        case IB_CACHED_FREE_FREE:
        {
            uint32_t vertex0FifoIndex = input.Decode( VertexDecoding, VERTEX_MAX_CODE_LENGTH );
            uint32_t relativeVertex1  = input.ReadVInt();
            uint32_t relativeVertex2  = input.ReadVInt();

            triangle[ 0 ]                                         = static_cast< IndiceType >( vertexFifo[ ( ( verticesRead - 1 ) - vertex0FifoIndex ) & VERTEX_FIFO_MASK ] );

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 1 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex1 );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex2 );

            verticesRead += 2;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;

            break;
        }
        case IB_FREE_FREE_FREE:
        {
            uint32_t relativeVertex0 = input.ReadVInt();
            uint32_t relativeVertex1 = input.ReadVInt();
            uint32_t relativeVertex2 = input.ReadVInt();

            vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         =
            triangle[ 0 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex0 );
            vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] =
            triangle[ 1 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex1 );
            vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] =
            triangle[ 2 ]                                         = static_cast< IndiceType >( ( newVertices - 1 ) - relativeVertex2 );

            verticesRead += 3;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

            ++edgesRead;

            break;
        }		
        }

        edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 1 ], triangle[ 2 ], triangle[ 0 ] );

        ++edgesRead;

        edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 2 ], triangle[ 0 ], triangle[ 1 ] );

        ++edgesRead;
    }

    // Skip over padding at the end, put in so a short prefix code (1 bit) doesn't cause an overflow
    input.Read( 32 );
}

// 32 bit indice/32bit attribute decompression
void DecompressMesh( 
    uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int32_t* vertexAttributes,
    ReadBitstream& input )
{
    DecompressMeshPrefix<uint32_t, int32_t>( triangles, triangleCount, vertexAttributeCount, vertexAttributes, input );
}

// 16 bit indice/32bit attribute decompression
void DecompressMesh(
    uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int32_t* vertexAttributes,
    ReadBitstream& input )
{
    DecompressMeshPrefix<uint16_t, int32_t>( triangles, triangleCount, vertexAttributeCount, vertexAttributes, input );
}

// 32 bit indice/32bit attribute decompression
void DecompressMesh(
    uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int16_t* vertexAttributes,
    ReadBitstream& input )
{
    DecompressMeshPrefix<uint32_t, int16_t>( triangles, triangleCount, vertexAttributeCount, vertexAttributes, input );
}

// 16 bit indice/32bit attribute decompression
void DecompressMesh(
    uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t vertexAttributeCount,
    int16_t* vertexAttributes,
    ReadBitstream& input )
{
    DecompressMeshPrefix<uint16_t, int16_t>( triangles, triangleCount, vertexAttributeCount, vertexAttributes, input );
}
