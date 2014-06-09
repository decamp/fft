/*
 * Copyright (c) 2014, Massachusetts Institute of Technology
 * Released under the BSD 2-Clause License
 * http://opensource.org/licenses/BSD-2-Clause
 */
package bits.fft;

import org.junit.Ignore;
import org.junit.Test;
import java.util.Random;


public class FastFourierTransformTest {

    private static final int OFFSET_0 = 3;
    private static final int DIM_0 = 16;

    private static final double[] INPUT_0 = {
            0.0,
            0.0,
            0.0,
            0.46193557475331404,
            -0.5189271686570283,
            0.2748348507002165,
            0.10087401023526787,
            0.19509055559440358,
            -0.33356320104670045,
            -0.229621630518563,
            0.9696830803996179,
            0.7583650357449603,
            0.8824983589642288,
            -0.45009206792903034,
            -0.7422056982524465,
            -0.7067966847069636,
            -0.9535237550322211,
            0.09347951439693114,
            0.9289737213537002,
            -0.7910186274980566,
            0.2502927269311186,
            -0.17840760901787656,
            0.552624582549865,
            0.981445571429566,
            -0.025534305939714397,
            0.49248281064466104,
            0.4663041403899877,
            0.6345941428186488,
            0.6777807000940366,
            0.05339886920973225,
            0.7986700232229871,
            -0.7321203188262155,
            -0.8338752035501702,
            0.9571486802956806,
            0.4447142383776974
    };

    private static final double[] OUTPUT_0 = {
            1.8147186670914086,
            2.6647862500402257,
            0.7682412719746914,
            0.7188424877600267,
            1.0782026349619975,
            -3.8579052912806806,
            2.4527663477593253,
            -2.31078432783205,
            -0.7734552604372538,
            5.051896414756831,
            2.0551633049643874,
            0.8590941482485186,
            -3.5898486789194743,
            0.616771471003068,
            -1.3023826623826356,
            0.05503965584308146,
            -0.21172816847209464,
            -4.374489946513126,
            0.3943882499269065,
            -1.6949732002136546,
            -1.6656841943867466,
            -5.0308279714019655,
            0.18676644876993675,
            2.362557260761443,
            3.4259692650934057,
            1.824385751045493,
            2.612894617620266,
            -3.454926589852724,
            -2.7108386868891836,
            0.9563077885428785,
            2.8557960393780872,
            -2.688608599419816
    };


    private static final int OFFSET_1 = 3;
    private static final int DIM_1 = 16;

    private static final double[] INPUT_1 = {
            0.0,
            0.0,
            0.0,
            0.4617563814065817,
            -0.17983837701559668,
            -0.5845703173805659,
            -0.33456588808097765,
            0.9355118188482414,
            -0.9877656354684774,
            0.9274095940464153,
            0.8797307775638197,
            0.8943898353263877,
            0.8741642977919393,
            -0.2056513156305888,
            -0.3049639415937795,
            -0.41188593599192647,
            0.012967254652470173,
            -0.7680658239346845,
            0.5410717601583555
    };

    private static final double[] OUTPUT_1 = {
            1.749694484697614,
            0.0,
            -2.814467725329596,
            -1.0526882912610054,
            3.4696444346603688,
            1.2261912793081786,
            -0.39597165158346614,
            0.6830510771354942,
            2.5106499624887086,
            1.0617451680870829,
            2.4643414658386105,
            -0.14985249308147097,
            -1.804603766907061,
            -0.6729395269375924,
            -0.9844359046047719,
            3.503999157882701,
            0.7480939886821072,
            0.0,
            -0.9844359046047715,
            -3.503999157882702,
            -1.8046037669070605,
            0.6729395269375924,
            2.4643414658386105,
            0.14985249308147164,
            2.5106499624887078,
            -1.0617451680870829,
            -0.3959716515834668,
            -0.6830510771354928,
            3.4696444346603705,
            -1.2261912793081786,
            -2.8144677253295964,
            1.0526882912610043
    };

    private static final double[] OUTPUT_1_INV = {
            0.10935590529360087,
            0.0,
            -0.17590423283309975,
            0.06579301820381284,
            0.21685277716627305,
            -0.07663695495676116,
            -0.024748228223966634,
            -0.042690692320968386,
            0.1569156226555443,
            -0.06635907300544268,
            0.15402134161491315,
            0.009365780817591936,
            -0.11278773543169131,
            0.042058720433599525,
            -0.061527244037798245,
            -0.2189999473676688,
            0.0467558742926317,
            0.0,
            -0.06152724403779822,
            0.21899994736766887,
            -0.11278773543169128,
            -0.042058720433599525,
            0.15402134161491315,
            -0.009365780817591977,
            0.15691562265554423,
            0.06635907300544268,
            -0.024748228223966676,
            0.0426906923209683,
            0.21685277716627316,
            0.07663695495676116,
            -0.17590423283309978,
            -0.06579301820381277
    };


    @Test
    @Ignore
    public void genMatrix() {
        int off = 3;
        int dim = 16;
        int len = dim;

        Random rand = new Random( 1 );
        double[] d = new double[off + len];

        for( int i = 0; i < len; i++ ) {
            d[i + off] = rand.nextDouble() * 2.0 - 1.0;
        }

        System.out.println( TestUtil.arrayToJava( d, 0, len + off, "INPUT_0" ) );
    }


    @Test
    public void testSpeed() {
        int len = 1024;

        double[] x = new double[len * 2];
        Random rand = new Random( 0 );

        for( int i = 0; i < len; i++ ) {
            x[i * 2] = rand.nextDouble() * 2.0 - 1.0;
            x[i * 2 + 1] = rand.nextDouble() * 2.0 - 1.0;
        }

        double[] out = new double[len * 2];
        FastFourierTransform trans = FastFourierTransform.create( len );

        Timer.start();

        for( int i = 0; i < 5000; i++ ) {
            trans.applyComplex( x, 0, false, out, 0 );
        }

        Timer.printSeconds( "FastFourierTransform time: " );
    }


    @Test
    public void testComplex() {
        final int dim = DIM_0;
        final int len = dim * 2;
        final int off = 7;

        double[] b = new double[len + off];
        double[] c = new double[len + off];

        FastFourierTransform trans = FastFourierTransform.create( dim );
        trans.applyComplex( INPUT_0, OFFSET_0, false, b, off );
        trans.applyComplex( b, off, true, c, off );

        //System.out.println(TestUtil.complexToMatlab(INPUT_0, OFFSET_0, DIM_0, 1));
        //System.out.println("===");
        //System.out.println(TestUtil.formatComplex(b, off, dim, 1));
        //System.out.println("===");
        //System.out.println(TestUtil.arrayToJava(b, off, len, "OUTPUT_0"));

        TestUtil.assertNear( OUTPUT_0, 0, b, off, len );
        TestUtil.assertNear( INPUT_0, OFFSET_0, c, off, len );
    }


    @Test
    public void testReal() {
        final int dim = DIM_1;
        final int len = dim * 2;
        final int off = 7;

        double[] b = new double[len + off];
        double[] c = new double[len + off];

        FastFourierTransform trans = FastFourierTransform.create( dim );
        trans.applyReal( INPUT_1, OFFSET_1, false, b, off );
        trans.applyReal( INPUT_1, OFFSET_1, true, c, off );

        //System.out.println(TestUtil.realToMatlab(INPUT_1, OFFSET_1, DIM_1, 1));
        //System.out.println("===");
        //System.out.println(TestUtil.formatComplex(c, off, dim, 1));
        //System.out.println("===");
        //System.out.println(TestUtil.arrayToJava(c, off, len, "OUTPUT_1_INV"));

        TestUtil.assertNear( OUTPUT_1, 0, b, off, len );
        TestUtil.assertNear( OUTPUT_1_INV, 0, c, off, len );
    }

}