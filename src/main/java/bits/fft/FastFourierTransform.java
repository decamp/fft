/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

/**
 * Performs a Fast Fourier Transform on an array of values.
 * Compatible with both real and complex valued inputs.
 * Outputs are always complex.
 * <p>
 * This class is reentrant (thread-safe).
 *
 * @author Philip DeCamp
 */
public class FastFourierTransform {

    private static final int MAX_BITS = 30;

    private final int mDim;
    private final int mBits;


    /**
     * The soled argument, <tt>dim</tt>, indicates the size
     * of vectors on which the transform will operate.  This
     * must be a power-of-two.
     * <p>
     * Memory footprint is about 18 bytes.
     *
     * @param dim Size of vector on which the transform operates.  Must be power-of-two.
     * @throws IllegalArgumentException if dim is smaller than 2, too large, or not a power-of-two.
     */
    public FastFourierTransform( int dim ) {
        mDim = dim;
        mBits = computeBitNum( dim );
    }


    /**
     * Performs a Fast Fourier Transform on a vector of complex values.
     * <p>
     * Complex samples must be stored in an array in a
     * tightly packed format: <br>
     * <tt>[... r0, i0, r1, i1, r2, i2 ...]</tt><br>
     * where <tt>rk</tt> is the real component of an lement and
     * <tt>ik</tt> is the corresponding imaginary component.
     *
     * @param x       Input array of complex samples: <b>NOTE:</b> <tt>x.length &gt= dim * 2 + xOff</tt>.
     * @param xOff    Start position of data in the input array.
     * @param inverse Set to false for normal FFT, true for inverse FFT.
     * @param out     Output array where transformed, complex elements are stored: <b>NOTE:</b> <tt>out.length &gt= dim * 2 + outOff</tt>
     * @param outOff  Start position into output array.
     */
    public void applyComplex( double[] x, int xOff, boolean inverse, double[] out, int outOff ) {
        // Bit-reverse the order of the data and copy into the output array.
        reverseBitShuffle( x, xOff, out, outOff, mDim, mBits );
        transform( out, outOff, mDim, inverse );

        //Perform scaling if this is an inverse transform.
        if( inverse ) {
            double scale = 1.0 / mDim;
            int len2 = mDim * 2;
            for( int i = 0; i < len2; i++ ) {
                out[i + outOff] *= scale;
            }
        }
    }

    /**
     * Performs a Fast Fourier Transform on an array of real values.
     * Note that the output samples are complex, so the output array
     * will need to hold twice as many double values as the input array.
     * <p>
     * Complex samples must be stored in an array in a
     * tightly packed format: <br>
     * <tt>[... r0, r1, r2, ...]</tt><br>
     * where <en>rk</en> is the real component of a sample.
     *
     * @param x          Input array of real-valued samples: <b>NOTE:</b> <tt>x.length &gt= dim + xOff</tt>.
     * @param xOff       Start position of data in the input array.
     * @param inverse    Set to false for normal FFT, true for inverse FFT.
     * @param out        Output array where transformed, complex samples are stored: <b>NOTE:</b> <tt>out.length &gt= numSamples * 2 + outOff</tt>
     * @param outOff     Start position into output array.
     */
    public void applyReal( double[] x, int xOff, boolean inverse, double[] out, int outOff ) {
        // Bit-reverse the order of the data and copy into the output array.
        final int shift = 31 - mBits;
        final int dim   = mDim;

        for( int i = 0; i < dim; i++ ) {
            int ii = i + xOff;
            int jj = (reverse( i ) >>> shift) + outOff;
            out[jj    ] = x[ii];
            out[jj + 1] = 0;
        }

        transform( out, outOff, mDim, inverse );

        // Perform scaling if this is an inverse transform.
        if( inverse ) {
            final double scale = 1.0 / mDim;
            final int len2 = mDim * 2;

            for( int i = 0; i < len2; i++ ) {
                out[i + outOff] *= scale;
            }
        }

    }



    static int computeBitNum( int n ) {
        int bits = 31 - Integer.numberOfLeadingZeros( n );
        if( bits <= 0 || bits >= MAX_BITS || 1 << bits != n ) {
            throw new IllegalArgumentException( "Dimension must be a power of two, larger than 1, and smaller than 1 << " + MAX_BITS );
        }
        return bits;
    }

    /**
     * @return <tt>val</tt> with bits in reversed order.
     */
    static int reverse( int val ) {
        long v = (((((val >>> 24)       ) * 0x0202020202L & 0x010884422010L) % 1023L)      ) |
                 (((((val >>> 16) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 8 ) |
                 (((((val >>>  8) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 16) |
                 (((((val       ) & 0xFF) * 0x0202020202L & 0x010884422010L) % 1023L) << 24);

        return (int)v;
    }



    static void transform( double[] x, int off, int len, boolean inverse ) {
        final double sign = inverse ? -1.0 : 1.0;

        int blockEnd = 1;
        double ar0, ar1, ar2;
        double ai0, ai1, ai2;

        for( int blockSize = 2; blockSize <= len; blockSize <<= 1 ) {
            final double angle = Math.PI * 2.0 / blockSize;
            final double cm1 = Math.cos( angle );

            //Use trig identity to quickly compute: sm1 = Math.sin(angle)
            //Block size >= 2  ergo   0 <= angle <= Math.PI  ergo  sin(angle) >= 0
            final double sm1 = sign * Math.sqrt( 1.0 - cm1 * cm1 );

            //Use trig identity to quickly compute: cm2 = Math.cos(2.0 * angle)
            final double cm2 = 2.0 * cm1 * cm1 - 1.0;

            //Use trig identity to quickly compute: sm2 = Math.sin(2.0 * angle);
            final double sm2 = 2.0 * sm1 * cm1;

            final double w = 2.0 * cm1;

            for( int i = 0; i < len; i += blockSize ) {
                ar2 = cm2;
                ar1 = cm1;

                ai2 = sm2;
                ai1 = sm1;

                for( int j = i * 2 + off, n = 0; n < blockEnd; j += 2, n++ ) {
                    ar0 = w * ar1 - ar2;
                    ar2 = ar1;
                    ar1 = ar0;

                    ai0 = w * ai1 - ai2;
                    ai2 = ai1;
                    ai1 = ai0;

                    int k = j + blockEnd * 2;
                    double tr = ar0 * x[k    ] - ai0 * x[k + 1];
                    double ti = ar0 * x[k + 1] + ai0 * x[k    ];

                    x[k    ] = x[j    ] - tr;
                    x[k + 1] = x[j + 1] - ti;

                    x[j    ] += tr;
                    x[j + 1] += ti;
                }
            }

            blockEnd = blockSize;
        }
    }


    private static void reverseBitShuffle( double[] a, int offA, double[] b, int offB, int dim, int bits ) {
        final int shift = 31 - bits;
        for( int xa = 0; xa < dim; xa++ ) {
            int xb = ( reverse( xa ) >>> shift ) + offB;
            b[xb    ] = a[offA++];
            b[xb + 1] = a[offA++];
        }
    }

}