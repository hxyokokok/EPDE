package hxy;

import hxy.DECOMPOSITION_METHOD;
import org.apache.commons.math3.exception.NumberIsTooLargeException;
import org.apache.commons.math3.random.RandomDataGenerator;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.pseudorandom.JMetalRandom;

import java.io.BufferedReader;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.util.StringTokenizer;

public class Tools {
    private static RandomDataGenerator rndGen;

    public static RandomDataGenerator getRndGen() {
        if (rndGen == null) {
            rndGen = new RandomDataGenerator() {
                @Override
                public void reSeed(long seed) {
                    reSeed(seed, seed);
                }

                public void reSeed(long JMetalSeed, long javaSeed) {
                    JMetalRandom.getInstance().setSeed(JMetalSeed);
                    super.reSeed(javaSeed);
                }

                @Override
                public int nextInt(int lower, int upper) throws NumberIsTooLargeException {
                    return lower == upper ? lower : super.nextInt(lower, upper);
                }
            };
            rndGen.reSeed(System.currentTimeMillis());
        }
        return rndGen;
    }

    public static void WS_transform(double lambda[][]) {
        for (int i = 0; i < lambda.length; i++) {
            double sum = 0;
            for (int j = 0; j < lambda[0].length; j++) {
                lambda[i][j] = lambda[i][j] > 1e-6 ? 1.0 / lambda[i][j] : 1e6;
                sum += lambda[i][j];
            }
            for (int j = 0; j < lambda[0].length; j++) {
                lambda[i][j] /= sum;
                if (lambda[i][j] < 1e-6) {
                    lambda[i][j] = 1e-6;
                }
            }
        }
    }

    public static double fitnessFunction(DoubleSolution individual, double[] idealPoint, double[] lambda, DECOMPOSITION_METHOD decompositionMethod) throws JMetalException {
        return fitnessFunction(individual, idealPoint, lambda, decompositionMethod, null);
    }

    public static double fitnessFunction(DoubleSolution individual, double[] idealPoint, double[] lambda, DECOMPOSITION_METHOD decompositionMethod, double[] z_nad) throws JMetalException {
        double fitness;
        switch (decompositionMethod) {
            case ASF:
                double maxFun = -1.0e+30;
                for (int n = 0; n < lambda.length; n++) {
                    double diff = individual.getObjective(n) - idealPoint[n];
                    if (z_nad != null) diff /= (z_nad[n] - idealPoint[n]);
                    double feval = diff / lambda[n];
                    if (feval > maxFun) {
                        maxFun = feval;
                    }
                }
                fitness = maxFun;
                break;
            case TCHE:
                maxFun = -1.0e+30;
                for (int n = 0; n < lambda.length; n++) {
                    double diff = Math.abs(individual.getObjective(n) - idealPoint[n]);
                    if (z_nad != null) {
                        diff /= (z_nad[n] - idealPoint[n]);
                    }
                    double feval;
                    if (lambda[n] == 0) {
                        feval = 0.0001 * diff;
                    } else {
                        feval = diff * lambda[n];
                    }
                    if (feval > maxFun) {
                        maxFun = feval;
                    }
                }
                fitness = maxFun;
                break;
            case WS:
                double sum = 0.0;
                for (int n = 0; n < lambda.length; n++) {
                    double diff = (individual.getObjective(n) - idealPoint[n]);
                    if (z_nad != null) {
                        diff /= (z_nad[n] - idealPoint[n]);
                    }
                    sum += (lambda[n]) * diff;
                }
                fitness = sum;
                break;
            case WP:
                sum = 1;
                for (int n = 0; n < lambda.length; n++) {
                    double diff = (individual.getObjective(n) - idealPoint[n]);
                    if (z_nad != null) {
                        diff /= (z_nad[n] - idealPoint[n]);
                    }
                    sum *= Math.pow(diff, lambda[n]);
                }
                fitness = sum;
                break;
            case PBI:
                double d1, d2, nl;
                double theta = 5.0;
                d1 = d2 = nl = 0.0;
                for (int i = 0; i < lambda.length; i++) {
                    double diff = (individual.getObjective(i) - idealPoint[i]);
                    if (z_nad != null) {
                        diff /= (z_nad[i] - idealPoint[i]);
                    }
                    d1 += diff * lambda[i];
                    nl += Math.pow(lambda[i], 2.0);
                }
                nl = Math.sqrt(nl);
                d1 = Math.abs(d1) / nl;
                for (int i = 0; i < lambda.length; i++) {
                    double diff = (individual.getObjective(i) - idealPoint[i]);
                    if (z_nad != null) {
                        diff /= (z_nad[i] - idealPoint[i]);
                    }
                    d2 += Math.pow(diff - d1 * (lambda[i] / nl), 2.0);
                }
                d2 = Math.sqrt(d2);
                fitness = (d1 + theta * d2);
                break;
            case SUM:
                fitness = 0;
                for (int i = 0; i < individual.getNumberOfObjectives(); i++) {
                    double diff = (individual.getObjective(i) - idealPoint[i]);
                    if (z_nad != null) {
                        diff /= (z_nad[i] - idealPoint[i]);
                    }
                    fitness += diff;
                }
                break;
            default:
                throw new JMetalException(" MOEAD.fitnessFunction: unknown type " + decompositionMethod);
        }
        return fitness;
    }

    public static double[][] ReadWeights(int M, int H) {
//        double lambda[][] = null;
//        switch (M) {
//            case 3:
//                lambda = new TwoLevelWeightVectorGenerator(12, 0, M).getLambda();
//                break;
//            case 5:
//                lambda = new TwoLevelWeightVectorGenerator(6, 0, M).getLambda();
//                break;
//            case 8:
//                lambda = new TwoLevelWeightVectorGenerator(3, 2, M).getLambda();
//                break;
//            case 10:
//                lambda = new TwoLevelWeightVectorGenerator(3, 2, M).getLambda();
//                break;
//            case 15:
//                lambda = new TwoLevelWeightVectorGenerator(2, 1, M).getLambda();
//                break;
//        }
        double lambda[][] = new double[H][M];
        if ((M == 2) && (H <= 300)) {
            for (int n = 0; n < H; n++) {
                double a = 1.0 * n / (H - 1);
                lambda[n][0] = a;
                lambda[n][1] = 1 - a;
            }
        } else {
            String dataFileName;
            dataFileName = M + "D_" + H + ".txt";
            try {
                InputStream in = Tools.class.getResourceAsStream("/weights_/" + dataFileName);
                InputStreamReader isr = new InputStreamReader(in);
                BufferedReader br = new BufferedReader(isr);
                int i = 0;
                int j = 0;
                String aux = br.readLine();
                while (aux != null) {
                    StringTokenizer st = new StringTokenizer(aux);
                    j = 0;
                    while (st.hasMoreTokens()) {
                        double value = new Double(st.nextToken());
                        lambda[i][j] = value;
                        j++;
                    }
                    aux = br.readLine();
                    i++;
                }
                br.close();
            } catch (Exception e) {
                throw new JMetalException("initializeUniformWeight: failed when reading file: ", e);
            }
        }
        for (int i = 0; i < lambda.length; i++) {
            for (int j = 0; j < lambda[0].length; j++) {
                if (lambda[i][j] < 1e-6) {
                    lambda[i][j] = 1e-6;
                }
            }
        }
        return lambda;
    }
}
