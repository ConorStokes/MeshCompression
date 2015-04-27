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
#include "meshcompression.h"
#include "writebitstream.h"
#include "indexcompressionconstants.h"
#include "meshcompressionconstants.h"
#include <assert.h>

#ifdef _MSC_VER
#define MC_INLINE __forceinline
#else
#define MC_INLINE inline
#endif 

#include "indexbufferencodetables.h"

// Classify a vertex as new, cached or free, outputting the relative position in the vertex indice cache FIFO.
static MC_INLINE VertexClassification ClassifyVertex( uint32_t vertex, const uint32_t* vertexRemap, const uint32_t* vertexFifo, uint32_t verticesRead, uint32_t& cachedVertexIndex )
{
    if ( vertexRemap[ vertex ] == VERTEX_NOT_MAPPED )
    {
        return NEW_VERTEX;
    }
    else
    {
        int32_t lowestVertexCursor = verticesRead >= VERTEX_FIFO_SIZE ? verticesRead - VERTEX_FIFO_SIZE : 0;

        // Probe backwards in the vertex FIFO for a cached vertex
        for ( int32_t vertexCursor = verticesRead - 1; vertexCursor >= lowestVertexCursor; --vertexCursor )
        {
            if ( vertexFifo[ vertexCursor & VERTEX_FIFO_MASK ] == vertex )
            {
                cachedVertexIndex = ( verticesRead - 1 ) - vertexCursor;

                return CACHED_VERTEX;
            }
        }

        return FREE_VERTEX;
    }
}


// Compress using triangle codes/prefix coding.
template <typename IndiceType, typename AttributeType>
void CompressMesh(
    const IndiceType* triangles,
    uint32_t triangleCount,
    uint32_t* vertexRemap,
    uint32_t vertexCount,
    uint32_t vertexAttributeCount,
    const AttributeType* vertexAttributes,
    WriteBitstream& output )
{
    EdgeTriangle edgeFifo[ EDGE_FIFO_SIZE ];
    uint32_t     vertexFifo[ VERTEX_FIFO_SIZE ];

    uint32_t          edgesRead    = 0;
    uint32_t          verticesRead = 0;
    uint32_t          newVertices  = 0;
    const IndiceType* triangleEnd  = triangles + ( triangleCount * 3 );

    assert( vertexCount < 0xFFFFFFFF );

    uint32_t* vertexRemapEnd = vertexRemap + vertexCount;

    // array of exponential moving average values to estimate optimal k for exp golomb codes
    // note we use 16/16 unsigned fixed point.
    uint32_t kArray[ 64 ];

    for ( uint32_t vertexAttributeIndex = 0; vertexAttributeIndex < vertexAttributeCount; ++vertexAttributeIndex )
    {
        // prime the array of ks for exp golomb with an average bitsize of 4
        // note that k is 
        kArray[ vertexAttributeIndex ] = 4 << 16;
    }

    // clear the vertex remapping to "not found" value of 0xFFFFFFFF - dirty, but low overhead.
    for ( uint32_t* remappedVertex = vertexRemap; remappedVertex < vertexRemapEnd; ++remappedVertex )
    {
        *remappedVertex = VERTEX_NOT_MAPPED;
    }

    // iterate through the triangles
    for ( const IndiceType* triangle = triangles; triangle < triangleEnd; triangle += 3 )
    {
        int32_t lowestEdgeCursor = edgesRead >= EDGE_FIFO_SIZE ? edgesRead - EDGE_FIFO_SIZE : 0;
        int32_t edgeCursor       = edgesRead - 1;
        bool    foundEdge        = false;

        int32_t spareVertex = 0;

        // check to make sure that there are no degenerate triangles.
        assert( triangle[ 0 ] != triangle[ 1 ] && triangle[ 1 ] != triangle[ 2 ] && triangle[ 2 ] != triangle[ 0 ] );

        // Probe back through the edge fifo to see if one of the triangle edges is in the FIFO
        for ( ; edgeCursor >= lowestEdgeCursor; --edgeCursor )
        {
            const EdgeTriangle& edge = edgeFifo[ edgeCursor & EDGE_FIFO_MASK ];

            // check all the edges in order and save the free vertex.
            if ( edge.second == triangle[ 0 ] && edge.first == triangle[ 1 ] )
            {
                foundEdge   = true;
                spareVertex = 2;
                break;
            }
            else if ( edge.second == triangle[ 1 ] && edge.first == triangle[ 2 ] )
            {
                foundEdge   = true;
                spareVertex = 0;
                break;
            }
            else if ( edge.second == triangle[ 2 ] && edge.first == triangle[ 0 ] )
            {
                foundEdge   = true;
                spareVertex = 1;
                break;
            }
        }

        // we found an edge so write it out, so classify a vertex and then write out the correct code.
        if ( foundEdge )
        {
            uint32_t cachedVertex;

            uint32_t             spareVertexIndice = triangle[ spareVertex ];
            VertexClassification freeVertexClass   = ClassifyVertex( spareVertexIndice, vertexRemap, vertexFifo, verticesRead, cachedVertex );
            uint32_t             relativeEdge      = ( edgesRead - 1 ) - edgeCursor;

            switch ( freeVertexClass )
            {
            case NEW_VERTEX:
            {
                output.WritePrefixCode( IB_EDGE_NEW, TrianglePrefixCodes );
                output.WritePrefixCode( relativeEdge, EdgePrefixCodes );

                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] = spareVertexIndice;
                vertexRemap[ spareVertexIndice ]              = newVertices;

                const EdgeTriangle&  edge               = edgeFifo[ edgeCursor & EDGE_FIFO_MASK ];
                const AttributeType* adjacent1Attribute = vertexAttributes + ( edge.first * vertexAttributeCount );
                const AttributeType* adjacent2Attribute = vertexAttributes + ( edge.second * vertexAttributeCount );
                const AttributeType* opposingAttribute  = vertexAttributes + ( edge.third * vertexAttributeCount );
                const AttributeType* vertexAttribute    = vertexAttributes + ( spareVertexIndice * vertexAttributeCount );
                const AttributeType* endAttributes      = vertexAttribute + vertexAttributeCount;
                uint32_t*            k                  = kArray;

                for ( ; vertexAttribute < endAttributes; ++adjacent1Attribute, ++adjacent2Attribute, ++opposingAttribute, ++vertexAttribute, ++k )
                {
                    int32_t  predicted = int32_t( *adjacent2Attribute ) + ( int32_t( *adjacent1Attribute ) - int32_t( *opposingAttribute ) );
                    int32_t  delta     = *vertexAttribute - predicted;
                    uint32_t kEstimate = output.WriteUniversalZigZag( delta, *k >> 16 );

                    // fixed point exponential moving average with alpha 0.125 (equivalent to N being 31)
                    *k = ( *k * 7 + ( kEstimate << 16 ) ) >> 3;
                }

                ++verticesRead;
                ++newVertices;
                break;
            }
            case CACHED_VERTEX:

                output.WritePrefixCode( IB_EDGE_CACHED, TrianglePrefixCodes );
                output.WritePrefixCode( relativeEdge, EdgePrefixCodes );
                output.WritePrefixCode( cachedVertex, CachedVertexPrefixCodes );

                break;

            case FREE_VERTEX:

                output.WritePrefixCode( IB_EDGE_FREE, TrianglePrefixCodes );
                output.WritePrefixCode( relativeEdge, EdgePrefixCodes );

                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] = spareVertexIndice;

                ++verticesRead;

                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ spareVertexIndice ] );

                break;
            }

            // Populate the edge fifo with the remaining edges
            // Note - the winding order is important as we'll need to re-produce this on decompression.
            // The edges are put in as if the found edge is the first edge in the triangle (which it will be when we
            // reconstruct).
            switch ( spareVertex )
            {
            case 0:

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 2 ], triangle[ 0 ], triangle[ 1 ] );

                ++edgesRead;

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

                ++edgesRead;
                break;

            case 1:

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 0 ], triangle[ 1 ], triangle[ 2 ] );

                ++edgesRead;

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 1 ], triangle[ 2 ], triangle[ 0 ] );

                ++edgesRead;
                break;

            case 2:

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 1 ], triangle[ 2 ], triangle[ 0 ] );

                ++edgesRead;

                edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( triangle[ 2 ], triangle[ 0 ], triangle[ 1 ] );

                ++edgesRead;
                break;
            }
        }
        else
        {
            VertexClassification classifications[ 3 ];
            uint32_t             cachedVertexIndices[ 3 ];

            // classify each vertex as new, cached or free, potentially extracting a cached indice.
            classifications[ 0 ] = ClassifyVertex( triangle[ 0 ], vertexRemap, vertexFifo, verticesRead, cachedVertexIndices[ 0 ] );
            classifications[ 1 ] = ClassifyVertex( triangle[ 1 ], vertexRemap, vertexFifo, verticesRead, cachedVertexIndices[ 1 ] );
            classifications[ 2 ] = ClassifyVertex( triangle[ 2 ], vertexRemap, vertexFifo, verticesRead, cachedVertexIndices[ 2 ] );

            // use the classifications to lookup the matching compression code and potentially rotate the order of the vertices.
            const VertexCompressionCase& compressionCase = CompressionCase[ classifications[ 0 ] ][ classifications[ 1 ] ][ classifications[ 2 ] ];

            // rotate the order of the vertices based on the compression classification.
            uint32_t reorderedTriangle[ 3 ];

            reorderedTriangle[ 0 ] = triangle[ compressionCase.vertexOrder[ 0 ] ];
            reorderedTriangle[ 1 ] = triangle[ compressionCase.vertexOrder[ 1 ] ];
            reorderedTriangle[ 2 ] = triangle[ compressionCase.vertexOrder[ 2 ] ];

            output.WritePrefixCode( compressionCase.code, TrianglePrefixCodes );

            switch ( compressionCase.code )
            {
            case IB_NEW_NEW_NEW:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = triangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = triangle[ 1 ];
                vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] = triangle[ 2 ];

                vertexRemap[ triangle[ 0 ] ] = newVertices;
                vertexRemap[ triangle[ 1 ] ] = newVertices + 1;
                vertexRemap[ triangle[ 2 ] ] = newVertices + 2;

                // interleaved verts
                // encode vert 0 absolute
                // encode vert 1 relative 0
                // encode vert 2 relative 0

                const AttributeType*  vert0    = vertexAttributes + ( triangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert1    = vertexAttributes + ( triangle[ 1 ] * vertexAttributeCount );
                const AttributeType*  vert2    = vertexAttributes + ( triangle[ 2 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++vert2, ++k )
                {
                    int32_t readVert0  = *vert0;

                    output.WriteUniversalZigZag( readVert0, EXP_GOLOMB_FIRST_NEW_K ); 

                    int32_t deltaVert1 = *vert1 - readVert0;

                    output.WriteUniversalZigZag( deltaVert1, ( *k >> 16 ) );

                    int32_t deltaVert2 = *vert2 - readVert0;

                    output.WriteUniversalZigZag( deltaVert2, ( *k >> 16 ) );
                }

                verticesRead += 3;
                newVertices  += 3;

                break;
            }

            case IB_NEW_NEW_CACHED:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 1 ];

                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 2 ] ], CachedVertexPrefixCodes );

                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;
                vertexRemap[ reorderedTriangle[ 1 ] ] = newVertices + 1;

                // encode vert 0 and 1 relative vert 2
                const AttributeType*  vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert1    = vertexAttributes + ( reorderedTriangle[ 1 ] * vertexAttributeCount );
                const AttributeType*  vert2    = vertexAttributes + ( reorderedTriangle[ 2 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++vert2, ++k )
                {
                    int32_t  readVert2 = *vert2;
                    int32_t deltaVert0 = *vert0 - readVert2;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );

                    int32_t deltaVert1 = *vert1 - readVert2;

                    output.WriteUniversalZigZag( deltaVert1, ( *k >> 16 ) );
                }

                verticesRead += 2;
                newVertices  += 2;

                break;
            }

            case IB_NEW_NEW_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 1 ];
                vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;
                vertexRemap[ reorderedTriangle[ 1 ] ] = newVertices + 1;

                // encode vert 0 and 1 relative vert 2
                const AttributeType*  vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert1    = vertexAttributes + ( reorderedTriangle[ 1 ] * vertexAttributeCount );
                const AttributeType*  vert2    = vertexAttributes + ( reorderedTriangle[ 2 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++vert2, ++k )
                {
                    int32_t readVert2  = *vert2;
                    int32_t deltaVert0 = *vert0 - readVert2;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );

                    int32_t deltaVert1 = *vert1 - readVert2;

                    output.WriteUniversalZigZag( deltaVert1, ( *k >> 16 ) );
                }

                verticesRead += 3;
                newVertices  += 2;

                break;
            }

            case IB_NEW_CACHED_CACHED:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] = reorderedTriangle[ 0 ];

                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 1 ] ], CachedVertexPrefixCodes );
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 2 ] ], CachedVertexPrefixCodes );
                
                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;

                // encode vert 0 relative vert 1
                const AttributeType*  vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert1    = vertexAttributes + ( reorderedTriangle[ 1 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;
                
                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++k )
                {
                    int32_t readVert1  = *vert1;
                    int32_t deltaVert0 = *vert0 - readVert1;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );
                }

                verticesRead += 1;
                newVertices  += 1;

                break;
            }
            case IB_NEW_CACHED_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 1 ] ], CachedVertexPrefixCodes );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;

                // encode vert 0 relative vert 1
                const AttributeType*  vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert1    = vertexAttributes + ( reorderedTriangle[ 1 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++k )
                {
                    int32_t readVert1  = *vert1;
                    int32_t deltaVert0 = *vert0 - readVert1;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );
                }

                verticesRead += 2;
                newVertices  += 1;

                break;
            }
            case IB_NEW_FREE_CACHED:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 1 ];

                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 1 ] ] );
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 2 ] ], CachedVertexPrefixCodes );

                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;

                // encode vert 0 relative vert 2
                const AttributeType*  vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType*  vert2    = vertexAttributes + ( reorderedTriangle[ 2 ] * vertexAttributeCount );
                const AttributeType*  vert0End = vert0 + vertexAttributeCount;
                const uint32_t*       k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert2, ++k )
                {
                    int32_t readVert2  = *vert2;
                    int32_t deltaVert0 = *vert0 - readVert2;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );
                }

                verticesRead += 2;
                newVertices  += 1;

                break;
            }
            case IB_NEW_FREE_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 1 ];
                vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 1 ] ] );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                vertexRemap[ reorderedTriangle[ 0 ] ] = newVertices;

                // encode vert 0 relative vert 1
                const AttributeType* vert0    = vertexAttributes + ( reorderedTriangle[ 0 ] * vertexAttributeCount );
                const AttributeType* vert1    = vertexAttributes + ( reorderedTriangle[ 1 ] * vertexAttributeCount );
                const AttributeType* vert0End = vert0 + vertexAttributeCount;
                const uint32_t*      k        = kArray;

                for ( ; vert0 < vert0End; ++vert0, ++vert1, ++k )
                {
                    int32_t readVert1  = *vert1;
                    int32_t deltaVert0 = *vert0 - readVert1;

                    output.WriteUniversalZigZag( deltaVert0, ( *k >> 16 ) );
                }

                verticesRead += 3;
                newVertices  += 1;

                break;
            }
            case IB_CACHED_CACHED_CACHED:
            {
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 0 ] ], CachedVertexPrefixCodes );
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 1 ] ], CachedVertexPrefixCodes );
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 2 ] ], CachedVertexPrefixCodes );

                break;
            }
            case IB_CACHED_CACHED_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 0 ] ], CachedVertexPrefixCodes );
                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 1 ] ], CachedVertexPrefixCodes );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                verticesRead += 1;

                break;
            }
            case IB_CACHED_FREE_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 1 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WritePrefixCode( cachedVertexIndices[ compressionCase.vertexOrder[ 0 ] ], CachedVertexPrefixCodes );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 1 ] ] );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                verticesRead += 2;

                break;
            }
            case IB_FREE_FREE_FREE:
            {
                vertexFifo[ verticesRead & VERTEX_FIFO_MASK ]         = reorderedTriangle[ 0 ];
                vertexFifo[ ( verticesRead + 1 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 1 ];
                vertexFifo[ ( verticesRead + 2 ) & VERTEX_FIFO_MASK ] = reorderedTriangle[ 2 ];

                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 0 ] ] );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 1 ] ] );
                output.WriteVInt( ( newVertices - 1 ) - vertexRemap[ reorderedTriangle[ 2 ] ] );

                verticesRead += 3;
                break;
            }

            default: // IB_EDGE_NEW, IB_EDGE_CACHED, IB_EDGE_0_NEW, IB_EDGE_1_NEW
                break;
            }

            // populate the edge fifo with the 3 most recent edges
            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( reorderedTriangle[ 0 ], reorderedTriangle[ 1 ], reorderedTriangle[ 2 ] );

            ++edgesRead;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( reorderedTriangle[ 1 ], reorderedTriangle[ 2 ], reorderedTriangle[ 0 ] );

            ++edgesRead;

            edgeFifo[ edgesRead & EDGE_FIFO_MASK ].set( reorderedTriangle[ 2 ], reorderedTriangle[ 0 ], reorderedTriangle[ 1 ] );

            ++edgesRead;
        }
    }

    // Pad out the buffer to make sure we don't overflow when trying to read the bits for the last prefix code table lookup.
    output.Write( 0, 32 );
}

void CompressMesh(
    const uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t* vertexRemap,
    uint32_t vertexCount,
    uint32_t vertexAttributeCount,
    const int32_t* vertexAttributes,
    WriteBitstream& output )
{
    CompressMesh< uint16_t, int32_t >( triangles, triangleCount, vertexRemap, vertexCount, vertexAttributeCount, vertexAttributes, output );
}

void CompressMesh(
    const uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t* vertexRemap,
    uint32_t vertexCount,
    uint32_t vertexAttributeCount,
    const int32_t* vertexAttributes,
    WriteBitstream& output )
{
    CompressMesh< uint32_t, int32_t >( triangles, triangleCount, vertexRemap, vertexCount, vertexAttributeCount, vertexAttributes, output );
}

void CompressMesh(
    const uint16_t* triangles,
    uint32_t triangleCount,
    uint32_t* vertexRemap,
    uint32_t vertexCount,
    uint32_t vertexAttributeCount,
    const int16_t* vertexAttributes,
    WriteBitstream& output )
{
    CompressMesh< uint16_t, int16_t >( triangles, triangleCount, vertexRemap, vertexCount, vertexAttributeCount, vertexAttributes, output );
}

void CompressMesh(
    const uint32_t* triangles,
    uint32_t triangleCount,
    uint32_t* vertexRemap,
    uint32_t vertexCount,
    uint32_t vertexAttributeCount,
    const int16_t* vertexAttributes,
    WriteBitstream& output )
{
    CompressMesh< uint32_t, int16_t >( triangles, triangleCount, vertexRemap, vertexCount, vertexAttributeCount, vertexAttributes, output );
}
