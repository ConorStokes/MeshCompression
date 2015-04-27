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
#ifndef WRITE_BIT_STREAM_H__
#define WRITE_BIT_STREAM_H__
#pragma once

#include <stdint.h>
#include <stdlib.h>
#include <memory.h>
#include <assert.h>

#ifdef _MSC_VER
#define WBS_INLINE __forceinline
#else
#define WBS_INLINE inline
#endif

#if defined (_MSC_VER)

#include <intrin.h>
#pragma intrinsic(_BitScanReverse)

#endif

// Used for prefix coding tables.
struct PrefixCode
{
    uint32_t code;
    uint32_t bitLength;
};

// Very simple bitstream for writing that will grow to accommodate written bits.
class WriteBitstream
{
public:

    // Construct the bit stream with an initial buffer capacity - should be a multiple of 8 and > 0
    WriteBitstream( size_t initialBufferCapacity = 16 )
    {
        m_bufferCursor =
        m_buffer       = new uint8_t[ initialBufferCapacity ];
        m_bufferEnd    = m_buffer + initialBufferCapacity;
        m_size         = 0;
        m_bitsLeft     = 64;
        m_bitBuffer    = 0;
    }

    ~WriteBitstream()
    {
        delete[] m_buffer;
    }

    // Size in bits.
    size_t Size() const { return m_size; }

    // Write a number of bits to the stream.
    void Write( uint32_t value, uint32_t bitCount );

    // Write a V int to the stream.
    void WriteVInt( uint32_t value );

    // Get the size in bytes
    size_t ByteSize() const { return ( m_size + 7 ) >> 3; }

    // Finish writing by flushing the buffer.
    void Finish();

    // Get the raw data for this buffer.
    const uint8_t* RawData() const { return m_buffer; }

    // Write a prefix code from the coding table to the stream.
    template <typename Ty>
    void WritePrefixCode( Ty input, const PrefixCode* codes );

    // Write exponential golomb like universal code, with a range 0..2^31 - 1, 
    // k is the fixed number of least significant bits encoded (bits above k will be variable encoded, with a unary code prefix describing their number)
    // k can be 0 to 31
    // Returns the k that would provide the smallest encoding for this number.
    uint32_t WriteUniversal( uint32_t value, uint32_t k );

    // Write exponential golomb like universal code, with a range -2^30..2^30 - 1, same as above but uses zig zag encoding
    // Returns the k that would provide the smallest encoding for this number.
    uint32_t WriteUniversalZigZag( int32_t value, uint32_t k );

    // Input can not be 0. Calculates floor( log2( input ) ) 
    static uint32_t Log2( uint32_t input );

    // Encode a signed 32-bit integer to an unsigned 32bit integer using zig-zag encoding.
    static uint32_t EncodeZigZag( int32_t input );

private:

    // If we need to grow the buffer.
    void GrowBuffer();
    
    // Not copyable
    WriteBitstream( const WriteBitstream& );

    // Not assignable
    WriteBitstream& operator=( const WriteBitstream& );

    uint64_t  m_bitBuffer;
    size_t    m_size;
    uint8_t*  m_buffer;
    uint8_t*  m_bufferCursor;
    uint8_t*  m_bufferEnd;
    uint32_t  m_bitsLeft;
};


WBS_INLINE uint32_t WriteBitstream::EncodeZigZag( int32_t input )
{
    return static_cast< uint32_t >( ( input << 1 ) ^ ( input >> 31 ) );
}


WBS_INLINE uint32_t WriteBitstream::WriteUniversalZigZag( int32_t input, uint32_t k )
{
    return WriteUniversal( EncodeZigZag( input ), k );
}


WBS_INLINE uint32_t WriteBitstream::Log2( uint32_t input )
{
#if defined ( _MSC_VER )

    unsigned long result;

    _BitScanReverse( &result, input );

    return static_cast< uint32_t >( result );

#elif defined( __GNUC__ ) || defined( __clang__ )

    return static_cast< uint32_t >( 31 - __builtin_clz( static_cast< unsigned int >( m_bitBuffer ) ) );

#else

    uint32_t result = ( input > 0xFFFF ) << 4;

    input >>= result;

    // This code should be a branchless count for the trailing 0 bits on pretty much every platform.
    uint32_t shift0 = ( input > 0xFF ) << 3;
    
    input >>= shift0;
    result |= shift0;

    uint32_t shift1 = ( input > 0xF ) << 2;

    input >>= shift1;
    result |= shift1;

    uint32_t shift2 = ( input > 0x3 ) << 1;

    input >>= shift2;
    result |= shift2;
    result |= input >> 1;

    return result;

#endif

}


WBS_INLINE uint32_t WriteBitstream::WriteUniversal( uint32_t value, uint32_t k )
{
    assert( k < 32 );
    assert( value < 0x80000000 );

    uint32_t bits = Log2( ( value << 1 ) | 1 );

    if ( bits <= k )
    {
        Write( 1, 1 );
        Write( value, k );
    }
    else
    {
        uint32_t bitsMinusK = bits - k;

        Write( uint32_t( 1 ) << bitsMinusK, bitsMinusK + 1 );
        Write( value & ~( uint32_t( 1 ) << ( bits - 1 ) ), bits - 1 );
    }

    return bits;
}


WBS_INLINE void WriteBitstream::Write( uint32_t value, uint32_t bitCount )
{
    m_bitBuffer |= ( static_cast<uint64_t>( value ) << ( 64 - m_bitsLeft ) ) & ( m_bitsLeft == 0 ? 0 : 0xFFFFFFFFFFFFFFFF );

    if ( bitCount > m_bitsLeft )
    {
        if ( m_bufferCursor > m_bufferEnd - 7 )
        {
            GrowBuffer();
        }

        m_bufferCursor[ 0 ] = m_bitBuffer & 0xFF;
        m_bufferCursor[ 1 ] = ( m_bitBuffer >> 8 ) & 0xFF;
        m_bufferCursor[ 2 ] = ( m_bitBuffer >> 16 ) & 0xFF;
        m_bufferCursor[ 3 ] = ( m_bitBuffer >> 24 ) & 0xFF;
        m_bufferCursor[ 4 ] = ( m_bitBuffer >> 32 ) & 0xFF;
        m_bufferCursor[ 5 ] = ( m_bitBuffer >> 40 ) & 0xFF;
        m_bufferCursor[ 6 ] = ( m_bitBuffer >> 48 ) & 0xFF;
        m_bufferCursor[ 7 ] = ( m_bitBuffer >> 56 ) & 0xFF;

        m_bufferCursor += 8;

        m_bitBuffer = value >> ( m_bitsLeft );
        m_bitsLeft  = 64 - ( bitCount - m_bitsLeft );
    }
    else
    {
        m_bitsLeft -= bitCount;
    }

    m_size += bitCount;
}


WBS_INLINE void WriteBitstream::WriteVInt( uint32_t value )
{
    do
    {
        uint32_t lower7 = value & 0x7F;

        value >>= 7;

        Write( lower7 | ( value > 0 ? 0x80 : 0 ), 8 );

    } while ( value > 0 );
}


inline void WriteBitstream::Finish()
{
    if ( m_bufferCursor > m_bufferEnd - 8 )
    {
        GrowBuffer();
    }

    m_bufferCursor[ 0 ] = m_bitBuffer & 0xFF;
    m_bufferCursor[ 1 ] = ( m_bitBuffer >> 8 ) & 0xFF;
    m_bufferCursor[ 2 ] = ( m_bitBuffer >> 16 ) & 0xFF;
    m_bufferCursor[ 3 ] = ( m_bitBuffer >> 24 ) & 0xFF;
    m_bufferCursor[ 4 ] = ( m_bitBuffer >> 32 ) & 0xFF;
    m_bufferCursor[ 5 ] = ( m_bitBuffer >> 40 ) & 0xFF;
    m_bufferCursor[ 6 ] = ( m_bitBuffer >> 48 ) & 0xFF;
    m_bufferCursor[ 7 ] = ( m_bitBuffer >> 56 ) & 0xFF;

    m_bufferCursor += 8;
}


WBS_INLINE void WriteBitstream::GrowBuffer()
{
    size_t    bufferSize     = m_bufferEnd - m_buffer;
    size_t    newBufferSize  = bufferSize * 2;
    size_t    bufferPosition = m_bufferCursor - m_buffer;
    uint8_t*  newBuffer      = new uint8_t[ newBufferSize ];

    ::memcpy( reinterpret_cast<void*>( newBuffer ), reinterpret_cast<void*>( m_buffer ), bufferSize );

    delete[] m_buffer;

    m_buffer       = newBuffer;
    m_bufferCursor = m_buffer + bufferPosition;
    m_bufferEnd    = m_buffer + newBufferSize;
}


template <typename Ty>
WBS_INLINE void WriteBitstream::WritePrefixCode( Ty input, const PrefixCode* codes )
{
    const PrefixCode& code = codes[ input ];

    Write( code.code, code.bitLength );
}


#endif // -- WRITE_BIT_STREAM_H__
