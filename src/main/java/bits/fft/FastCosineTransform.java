/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

/**
 * Performs a Fast Discrete Cosine Transform on an array of real values.
 * Currently relies on FastFourierTransform, so it could be faster
 * if it used a more specialized FCT algorithm.
 * <p/>
 * Not thread safe. For parallel operation, create multiple FastCosineTransform
 * objects.
 *
 * @author Philip DeCamp
 */
public class FastCosineTransform {

    /**
     * Factory method to create FastCosineTransform objects.
     * The sole argument, <code>dim</code> indicates the size
     * of vectors on which the transform will operate.  This
     * must be a power-of-two.
     * <p/>
     * Memory allocation is just above dim*6*8 bytes.
     *
     * @param dim Size of vector on which this transform operates.  Must be power-of-two.
     * @return new FastCosineTransform instance.
     * @throws IllegalArgumentException if dim is smaller than 2, too large, or not a power-of-two.
     */
    public static FastCosineTransform create( int dim ) {
        return new FastCosineTransform( dim );
    }

    /**
     * Performs a Fast Discrete Cosine Transform on an array of real values.
     * <p/>
     * Not thread safe.
     *
     * @param a       Input array of real values with size [dim].
     * @param aOff    Offset into array <code>a</code>
     * @param inverse Set to <code>true</code> to perform inverse transform.
     * @param out     Output matrix where DCT coeffs are stored.  Must have space for <code>dim</code> values.
     * @param outOff  Offset into array <code>out</code>
     */
    public void apply( double[] a, int aOff, boolean inverse, double[] out, int outOff ) {
        if( !inverse ) {
            shuffle1( a, aOff, mDim, mBits, mWorkA );
            FastFourierTransform.doFft1d( mWorkA, 0, mDim, false );
            shuffle2( mWorkA, mWeight, mDim, out, outOff );
        } else {
            invShuffle1( a, aOff, mInvWeight, mDim, mBits, mWorkA );
            FastFourierTransform.doFft1d( mWorkA, 0, mDim, true );
            invShuffle2( mWorkA, mDim, mBits, out, outOff );
        }
    }



    private static final int MAX_BITS = 30;


    private final int mDim;
    private final int mBits;

    private final double[] mWeight;
    private final double[] mInvWeight;
    private final double[] mWorkA;


    private FastCosineTransform( int dim ) {
        mDim = dim;
        mBits = computeNumBits( dim );

        mWeight = new double[dim * 2];
        mInvWeight = new double[dim * 2];
        mWorkA = new double[dim * 2];

        computeWeightVector( dim, mWeight, 0 );
        invComputeWeightVector( dim, mInvWeight, 0 );
    }


    private static int computeNumBits( int n ) {
        int bits = 1; //Start at 1 because 2^0 is invalid.

        for(; bits < MAX_BITS; bits++ ) {
            if( (1 << bits) == n ) {
                break;
            }
        }

        if( bits == MAX_BITS ) {
            throw new IllegalArgumentException( "Dimension must be a power of two, larger than 1, and smaller than 1 << " + MAX_BITS );
        }

        return bits;
    }


    private static void computeWeightVector( int dim, double[] w, int offW ) {
        w[0] = 1.0;
        w[1] = 0.0;

        for( int i = 1; i < dim; i++ ) {
            double cos = Math.cos( -i * Math.PI * 0.5 / dim );
            double sin = -Math.sqrt( 1.0 - cos * cos );

            w[i * 2] = 2.0 * cos;
            w[i * 2 + 1] = 2.0 * sin;
        }
    }


    private static void invComputeWeightVector( int dim, double[] w, int offW ) {
        final double s = 1.0 / dim;
        w[0] = 1.0 * s;
        w[1] = 0.0;

        for( int i = 1; i < dim; i++ ) {
            double cos = Math.cos( i * Math.PI * 0.5 / dim );
            double sin = Math.sqrt( 1.0 - cos * cos );

            w[i * 2] = cos * s;
            w[i * 2 + 1] = sin * s;
        }
    }

    /**
     * 1. Drop complex components <br/>
     * 2. Butterfly shuffle <br/>
     * 3. Reverse-bit shuffle <br/>
     */
    private static void shuffle1( double[] a, int offA, int dim, int bits, double[] out ) {
        final int shift = 30 - bits;
        final int dim2 = dim * 2;

        for( int rowOut2 = 0; rowOut2 < dim2; rowOut2 += 2 ) {
            //Bit reversal shuffle.
            int rowIn = ( reverse( rowOut2 ) >>> shift );

            //Butterfly shuffle.
            rowIn = rowIn + (rowIn / dim) * (dim2 - 2 * rowIn - 1);

            out[rowOut2] = a[rowIn + offA];
            out[rowOut2 + 1] = 0.0;
        }
    }

    /**
     * 1. Apply weight vector.
     * 2. Drop imaginary components.
     */
    private static void shuffle2( double[] a, double[] w, int dim, double[] out, int offOut ) {
        final int dim2 = 2 * dim;

        for( int y2 = 0; y2 < dim2; y2 += 2 ) {
            out[offOut++] = a[y2] * w[y2] - a[y2 + 1] * w[y2 + 1];
        }
    }

    /**
     * 1. Apply weights, <br/>
     * 2. Reverse-bit shuffle
     */
    private static void invShuffle1( double[] a, int offA, double[] w, int dim, int bits, double[] out ) {
        final int shift = 31 - bits;

        for( int rowIn = 0; rowIn < dim; rowIn++ ) {
            int rowOut2 = ( reverse( rowIn ) >>> shift );
            double v = a[offA + rowIn];

            out[rowOut2] = v * w[rowIn * 2];
            out[rowOut2 + 1] = v * w[rowIn * 2 + 1];
        }
    }

    /**
     * 1. drop imaginary component<br/>
     * 2. reverse butterfly shuffle<br/>
     */
    private static void invShuffle2( double[] a, int dim, int bits, double[] out, int offOut ) {
        final int dim2 = dim * 2;

        for( int i2 = 0; i2 < dim2; i2 += 2 ) {
            //Reverse butterfly shuffle.
            final int rowO = i2 + (i2 / dim) * (dim2 - 1 - 2 * i2);
            out[rowO + offOut] = a[i2];
        }
    }


    /**
     * @return Version of val with reversed bits.
     */
    private static int reverse( int val ) {
        long v = (((((val >>> 24)       ) * 0x0202020202L & 0x010884422010L) % 1023L)      ) |
                 (((((val >>> 16) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 8 ) |
                 (((((val >>>  8) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 16) |
                 (((((val       ) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 24);

        return (int)v;
    }

}