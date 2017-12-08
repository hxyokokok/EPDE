package hxy;

import org.uma.jmetal.algorithm.multiobjective.moead.MOEAD;
import org.uma.jmetal.algorithm.multiobjective.moead.util.MOEADUtils;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.JMetalException;
import org.uma.jmetal.util.SolutionListUtils;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * Created by hxy on 2017/3/26.
 */
public class MOEAD_EP extends MOEAD {
    DECOMPOSITION_METHOD decompositionMethod = null;
    // evolutionary path and the updating rate
    double[] EPrate = null;
    ArrayList<ArrayList<DoubleSolution>> EP = new ArrayList<>();

    public MOEAD_EP(MOP p) {
        this(p, p.getH(), p.getMaxFEs(), DECOMPOSITION_METHOD.PBI);
        HistoryLength = neighborSize;
        EPrate = new double[4];
        for (int i = 0; i < 4; i++)
            EPrate[i] = 1 - 1.0 * i / 4;
    }

//    public static MOEAD_EP createInstance(MOP p) {
//        DECOMPOSITION_METHOD decompositionMethod = DECOMPOSITION_METHOD.PBI;
//        MOEAD_EP m = new MOEAD_EP(p, p.getH(), p.getMaxFEs(), decompositionMethod);
//        m.HistoryLength = m.neighborSize;
//        m.EPrate = new double[4];
//        for (int i = 0; i < 4; i++)
//            m.EPrate[i] = 1 - 1.0 * i / 4;
//        return m;
//    }

//    public static MOEAD_EP createInstance(MOP p, double historyLengthRate, int EPrateLength) {
//        DECOMPOSITION_METHOD decompositionMethod = DECOMPOSITION_METHOD.PBI;
//        MOEAD_EP m = new MOEAD_EP(p, p.getH(), p.getMaxFEs(), decompositionMethod);
//        m.HistoryLength = (int) (m.neighborSize * historyLengthRate);
//        m.EPrate = new double[EPrateLength];
//        for (int i = 0; i < EPrateLength; i++)
//            m.EPrate[i] = 1 - 1.0 * i / EPrateLength;
//        return m;
//    }

    private MOEAD_EP(Problem<DoubleSolution> problem, int N, int maxFEs, DECOMPOSITION_METHOD decompositionMethod) {
        super(problem, N, N, maxFEs, null, null, null, null, 0.9, Integer.MAX_VALUE, N / 10);
        initializePopulation();
        lambda = Tools.ReadWeights(problem.getNumberOfObjectives(), populationSize);
//        initializeUniformWeight();
        initializeNeighborhood();
        initializeIdealPoint();
        this.decompositionMethod = decompositionMethod;
    }

    @Override
    public List<DoubleSolution> getResult() {
        return SolutionListUtils.getNondominatedSolutions(population);
    }

    //    initialize evolution path
    private void initEP() {
        for (int i = 0; i < EPrate.length; i++) {
            ArrayList<DoubleSolution> tmp = new ArrayList<>();
            this.population.forEach(s -> tmp.add((DoubleSolution) s.copy()));
            EP.add(tmp);
        }
    }

    // parameters for DE
    double[][] HistoryF, HistoryCr;
    int HistoryLength = -1;
    int[] HistoryPtr;

    public static double[] generateParametersForDE(double[] H_F, double[] H_Cr) {
        double mean_f, mean_cr;
        double f, cr;

        int rndIdx = Tools.getRndGen().nextInt(0, H_F.length - 1);
        mean_f = H_F[rndIdx];
        mean_cr = H_Cr[rndIdx];
        do {
            f = mean_f + Tools.getRndGen().nextGaussian(0, 1) * 0.2;
        } while (f > 1 || f < 0);
        do {
            cr = mean_cr + Tools.getRndGen().nextGaussian(0, 1) * 0.1;
        } while (cr > 1 || cr < 0);
        double[] ret = new double[2];
        ret[0] = f;
        ret[1] = cr;
        return ret;
    }

    private void initPara() {
        HistoryF = new double[populationSize][HistoryLength];
        HistoryCr = new double[populationSize][HistoryLength];
        HistoryPtr = new int[populationSize];
        for (int i = 0; i < populationSize; i++) {
            Arrays.fill(HistoryF[i], 0.5);
            Arrays.fill(HistoryCr[i], 0.5);
        }
        Arrays.fill(HistoryPtr, 0);
    }

    public static DoubleSolution DE_operation(DoubleSolution future, DoubleSolution current, DoubleSolution r1, DoubleSolution r2, double F, double Cr) {
        DoubleSolution V = (DoubleSolution) current.copy();
        int jrand = Tools.getRndGen().nextInt(0, V.getNumberOfVariables() - 1);
        for (int j = 0; j < V.getNumberOfVariables(); j++) {
            if (Tools.getRndGen().nextUniform(0, 1) < Cr || j == jrand) {
                double y = future.getVariableValue(j) + F * (r1.getVariableValue(j) - r2.getVariableValue(j));
                if (y > current.getUpperBound(j)) {
//                    double r = hxy.Tools.getRndGen().nextUniform(0,1);
//                    y = r * current.getVariableValue(j) + (1-r) * current.getUpperBound(j);
                    y = current.getUpperBound(j);
                } else if (y < current.getLowerBound(j)) {
//                    double r = hxy.Tools.getRndGen().nextUniform(0,1);
//                    y = r * current.getVariableValue(j) + (1-r) * current.getLowerBound(j);
                    y = current.getLowerBound(j);
                }
                V.setVariableValue(j, y);
            }
        }
        return V;
    }

    @Override
    public void run() {
        PolynomialMutation PMOpt = new PolynomialMutation(1.0 / problem.getNumberOfVariables(), 20);
        initEP(); // initialize EP
        initPara(); // initialize parameters for DE 
        evaluations = populationSize;
        do {
            int[] permutation = new int[populationSize];
            MOEADUtils.randomPermutation(permutation, populationSize);

            for (int i = 0; i < populationSize; i++) {
                int subProblemId = permutation[i];
                NeighborType neighborType = chooseNeighborType();

//                select two neighbors to generate difference vector
                List<DoubleSolution> parents = parentSelection(subProblemId, neighborType);

//                generate a future solution of this subproblem
                double[] tmp = generateParametersForDE(HistoryF[subProblemId], HistoryCr[subProblemId]);
                // double[] tmp = generateParametersForDE(subProblemId);
                double F = tmp[0], Cr = tmp[1];
                DoubleSolution futureSolution = generateFutureSolution(EP, subProblemId, EPrate, (1.0 - F) * (1 - 1.0 * evaluations / maxEvaluations) + 1.01);
                DoubleSolution V = DE_operation(futureSolution, population.get(subProblemId), parents.get(0), parents.get(1), F, Cr);
//                PM
                DoubleSolution child = PMOpt.execute(V);

                problem.evaluate(child);
                evaluations++;

                updateIdealPoint(child);
                List<Integer> changedSubproblemList = replaceNeighborhood(child, subProblemId, neighborType);

                for (int spIdx : changedSubproblemList) {
//                    updateNeighborhood evolutionary path
                    updateEvolutionaryPath(spIdx);
//                    updateNeighborhood parameters for DE use success history
                    HistoryF[spIdx][HistoryPtr[spIdx]] = F;
                    HistoryCr[spIdx][HistoryPtr[spIdx]] = Cr;
                    HistoryPtr[spIdx]++;
                    if (HistoryPtr[spIdx] >= HistoryLength) {
                        HistoryPtr[spIdx] = 0;
                    }
                }
            }
        } while (evaluations < maxEvaluations);
    }

    List<Integer> replaceNeighborhood(DoubleSolution individual, int subProblemId, NeighborType neighborType) throws JMetalException {
        int size;
        int time;

        time = 0;
        List<Integer> changedSubproblemList = new ArrayList<>();
        if (neighborType == NeighborType.NEIGHBOR) {
            size = neighborhood[subProblemId].length;
        } else {
            size = population.size();
        }
        int[] perm = new int[size];

        MOEADUtils.randomPermutation(perm, size);

        for (int i = 0; i < size; i++) {
            int k;
            if (neighborType == NeighborType.NEIGHBOR) {
                k = neighborhood[subProblemId][perm[i]];
            } else {
                k = perm[i];
            }
            double neighbor_fit_, u_fit_;

            neighbor_fit_ = Tools.fitnessFunction(population.get(k), idealPoint, lambda[k], decompositionMethod);
            u_fit_ = Tools.fitnessFunction(individual, idealPoint, lambda[k], decompositionMethod);

            if (u_fit_ < neighbor_fit_) {
                population.set(k, (DoubleSolution) individual.copy());
                time++;
                changedSubproblemList.add(k);
            }
            if (time >= maximumNumberOfReplacedSolutions) {
                break;
            }
        }
        return changedSubproblemList;
    }

    public static double[] getCoefficient(double[] EPrate, double t) {
        double[] coeff = new double[EPrate.length];
        Arrays.fill(coeff, 1);
        for (int i = 0; i < EPrate.length; i++) {
            for (int j = 0; j < EPrate.length; j++) {
                if (i == j) continue;
                coeff[i] *= (t - EPrate[j]) / (EPrate[i] - EPrate[j]);
            }
        }
        return coeff;
    }

    public static DoubleSolution generateFutureSolution(ArrayList<ArrayList<DoubleSolution>> EP, int subProblemId, double[] EPrate, double time) {
        DoubleSolution s = (DoubleSolution) EP.get(0).get(0).copy();
        double[] coff = getCoefficient(EPrate, time);
        for (int j = 0; j < s.getNumberOfVariables(); j++) {
            double V = 0;
            for (int i = 0; i < EPrate.length; i++) {
                V += coff[i] * EP.get(i).get(subProblemId).getVariableValue(j);
                // if (V > s.getUpperBound(j)) V = s.getUpperBound(j);
                // if (V < s.getLowerBound(j)) V = s.getLowerBound(j);
            }
            s.setVariableValue(j, V);
        }
        return s;
    }

    //                updateNeighborhood EP
    void updateEvolutionaryPath(int subproblemId) {
        for (int j = 0; j < EPrate.length; j++) {
            for (int t = 0; t < problem.getNumberOfVariables(); t++) {
                double originalX = EP.get(j).get(subproblemId).getVariableValue(t);
                double newX = originalX * (1 - EPrate[j]) + population.get(subproblemId).getVariableValue(t) * EPrate[j];
                EP.get(j).get(subproblemId).setVariableValue(t, newX);
            }
        }
    }

    @Override
    public String getName() {
        String name = "MOEA/D-EP";
        return name;
    }

}

