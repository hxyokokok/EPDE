package hxy;

import org.uma.jmetal.problem.DoubleProblem;
import org.uma.jmetal.problem.impl.AbstractDoubleProblem;
import org.uma.jmetal.solution.DoubleSolution;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by hxy on 2017/3/25.
 */
public class MOP extends AbstractDoubleProblem {
    public static List<MOP> createManyObjectiveBenchmark(String funName) throws Exception {
        List<MOP> plist = new ArrayList<>();
        int[] MList = new int[]{3, 5, 8, 10, 15};
        if (funName.equals("WFG")) {
            for (int pNo = 1; pNo <= 9; pNo++) {
                for (int M : MList) {
                    plist.add(new MOP("WFG", pNo, M));
                }
            }
        } else if (funName.equals("DTLZ1-4")) {
            for (int pNo = 1; pNo <= 4; pNo++) {
                for (int M : MList) {
                    plist.add(new MOP("DTLZ", pNo, M));
                }
            }
        } else if (funName.equals("DTLZ5-7")) {
            for (int pNo = 5; pNo <= 7; pNo++) {
                for (int M : MList) {
                    plist.add(new MOP("DTLZ", pNo, M));
                }
            }
        } else if (funName.equals("MINUS_DTLZ1-4")) {
            for (int pNo = 1; pNo <= 4; pNo++) {
                for (int M : MList) {
                    plist.add(new MOP("MINUS_DTLZ", pNo, M));
                }
            }
        }
        return plist;
    }

    private DoubleProblem innerProblem;
//    int maxFEs, dim, H;
    private int maxFEs, dim, H;
    private int K, L;
    private String pfName;


    public MOP(String funName, int funNo, int M) throws Exception {
        switch (M) {
            case 3:
                H = 91;
                break;
            case 5:
                H = 210;
                break;
            case 8:
                H = 156;
                break;
            case 10:
                H = 275;
                break;
            case 15:
                H = 135;
                break;
        }
        if (funName.equals("WFG")) {
            switch (M) {
                case 3:
                    maxFEs = 36400;
                    break;
                case 5:
                    maxFEs = 157500;
                    break;
                case 8:
                    maxFEs = 234000;
                    break;
                case 10:
                    maxFEs = 550000;
                    break;
                case 15:
                    maxFEs = 405000;
                    break;
                default:
                    throw new Exception("No problem configuration found for M=" + M);
            }
            L = 20;
            K = 2 * (M - 1);
            dim = K + L;
            innerProblem = (DoubleProblem) Class.forName("org.uma.jmetal.problem.multiobjective.wfg.WFG" + funNo).getDeclaredConstructor(Integer.class, Integer.class, Integer.class).newInstance(K, L, M);
            pfName = "WFG" + funNo + "_" + M + ".txt";
        } else if (funName.equals("DTLZ")) {
            if (funNo == 1) K = 5;
            else if (funNo <= 6) K = 10;
            else if (funNo == 7) K = 20;
            else throw new Exception("No such problem with fno=" + funNo);
            int[] maxGen = null;
            switch (M) {
                case 3:
                    maxGen = new int[]{400, 250, 1000, 600, 250, 250, 1000};
                    break;
                case 5:
                    maxGen = new int[]{600, 350, 1000, 1000, 350, 350, 1000};
                    break;
                case 8:
                    maxGen = new int[]{750, 500, 1000, 1250, 500, 500, 1000};
                    break;
                case 10:
                    maxGen = new int[]{1000, 750, 1500, 2000, 750, 750, 1500};
                    break;
                case 15:
                    maxGen = new int[]{1500, 1000, 2000, 3000, 1000, 1000, 2000};
                    break;
                default:
                    throw new Exception("No problem configuration found for M=" + M);
            }
            maxFEs = maxGen[funNo - 1] * H;
            dim = K + M - 1;
            innerProblem = (DoubleProblem) Class.forName("org.uma.jmetal.problem.multiobjective.dtlz.DTLZ" + funNo).getDeclaredConstructor(Integer.class, Integer.class).newInstance(dim, M);
            pfName = "DTLZ" + funNo + "_" + M + ".txt";
        } else if (funName.startsWith("MINUS")) {
            String ff = funName.substring(6);
            MOP p = new MOP(ff, funNo, M) {
                @Override
                public void evaluate(DoubleSolution solution) {
                    super.evaluate(solution);
                    for (int j = 0; j < M; j++) {
                        solution.setObjective(j, solution.getObjective(j) * -1);
                    }
                }
            };
            innerProblem = p;
            pfName = "MINUS_" + p.pfName;
            dim = p.dim;
            H = p.H;
            maxFEs = p.maxFEs;
            L = p.L;
            K = p.K;
        } else {
            System.err.printf("no such benchmark test function:\t" + funName);
        }
        this.setName(funName + funNo + "_" + M + "M_" + maxFEs + "fes");
        this.setNumberOfObjectives(M);
        this.setNumberOfVariables(dim);
    }

    @Override
    public void evaluate(DoubleSolution solution) {
        innerProblem.evaluate(solution);
    }

    public int getMaxFEs() {
        return maxFEs;
    }

    public int getH() {
        return H;
    }

    public String getPfName() {
        return pfName;
    }

    @Override
    public Double getUpperBound(int index) {
        return innerProblem.getUpperBound(index);
    }

    @Override
    public Double getLowerBound(int index) {
        return innerProblem.getLowerBound(index);
    }

    public static void main(String args[]) throws Exception {
        MOP mop = new MOP("WFG", 3, 8);
        DoubleSolution s = mop.createSolution();
        double[] x = {
                0.0860364757097445,
                1.90247193451748,
                1.74249830920364,
                7.74792823292554,
                6.04722667860648,
                5.55020092848946,
                10.4542126262948,
                9.52661233845861,
                6.84944040224494,
                18.1565291419875,
                8.95621053035107,
                10.0921496049842,
                2.58961045589386,
                5.73976906231547,
                27.7492013639356,
                10.5548710243126,
                18.6191429738934,
                27.7466211319516,
                29.8027571665599,
                21.2930542571824,
                41.9892559406942,
                20.2878993972932,
                31.7911578432091,
                17.9556421796001,
                21.2565900479344,
                30.5776021877456,
                36.8607736897153,
                4.9922078073677,
                25.5698217738811,
                3.06593366321319,
                43.9818108924038,
                36.7036581083647,
                20.0850819751629,
                40.1152385903161};
        for (int i = 0; i < s.getNumberOfVariables(); i++) {
            s.setVariableValue(i, x[i]);
        }
        mop.evaluate(s);
        for (int i = 0; i < s.getNumberOfObjectives(); i++) {
            System.out.println(s.getObjective(i));
        }
    }
}
