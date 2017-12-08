package hxy;

import org.apache.commons.math3.linear.*;
import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.operator.MutationOperator;
import org.uma.jmetal.operator.impl.mutation.PolynomialMutation;
import org.uma.jmetal.problem.Problem;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.solution.impl.DefaultDoubleSolution;
import org.uma.jmetal.util.SolutionListUtils;
import org.uma.jmetal.util.solutionattribute.Ranking;
import org.uma.jmetal.util.solutionattribute.impl.DominanceRanking;

import java.util.*;
import java.util.stream.Collectors;

class ReferencePoint {
    double[] lambda;
    List<DoubleSolution> candidates = new LinkedList<>();
    List<DoubleSolution> members = new LinkedList<>();
    List<DoubleSolution> evolutionPath = new ArrayList<>();
    double[] HistoryF, HistoryCr;
    int HistoryPtr = 0;
}

@SuppressWarnings("serial")
public class NSGA3_EP implements Algorithm<List<DoubleSolution>> {
    private List<DoubleSolution> population;
    private Problem<DoubleSolution> problem;
    private int FEs, maxFEs;
    private int N, M, dim;
    private List<ReferencePoint> referencePoints = new Vector<>();
    private MutationOperator<DoubleSolution> mutationOperator;
    private double[] z_star, z_nad;
    private int T;
    private double[] EPrate = new double[]{1, 0.75, 0.5, 0.25};

    private int historyLength;

    public NSGA3_EP(MOP p) {
        this(p, 1.0);
    }

    public NSGA3_EP(MOP p, double HLR) {
        this.problem = p;
        this.maxFEs = p.getMaxFEs();
        this.mutationOperator = new PolynomialMutation(1.0 / p.getNumberOfVariables(), 20);
        this.M = p.getNumberOfObjectives();
        this.dim = p.getNumberOfVariables();
        this.N = p.getH();
        while (N % 4 > 0) N++;
        double[][] Lambda = Tools.ReadWeights(M, p.getH());
        for (double[] lambda : Lambda) {
            ReferencePoint rp = new ReferencePoint();
            rp.lambda = lambda;
            referencePoints.add(rp);
        }
//        init population
        population = new ArrayList<>();
        for (int i = 0; i < N; i++) {
            population.add(problem.createSolution());
        }
        population.forEach(problem::evaluate);
        FEs = N;
//        init EP
        int cnt = 0;
        for (ReferencePoint w : referencePoints) {
            for (int i = 0; i < EPrate.length; i++) {
                w.evolutionPath.add((DoubleSolution) population.get(cnt).copy());
            }
            cnt++;
        }
//        init association relationship
        cnt = 0;
        for (DoubleSolution s : population) {
            referencePoints.get(cnt).members.add(s);
            cnt = ++cnt % p.getH();
        }
//        normalize reference point
        for (ReferencePoint w : referencePoints) {
            double norm = 0;
            for (double v : w.lambda) norm += v * v;
            norm = Math.sqrt(norm);
            for (int i = 0; i < w.lambda.length; i++) w.lambda[i] /= norm;
        }
//        init parameter history
        this.historyLength = (int) (N * HLR / 10);
        for (ReferencePoint w : referencePoints) {
            w.HistoryF = new double[historyLength];
            w.HistoryCr = new double[historyLength];
            Arrays.fill(w.HistoryF, 0.5);
            Arrays.fill(w.HistoryCr, 0.5);
        }
//        init parameters for the initial population
//        population.forEach(s -> s.setAttribute("para", new double[]{0.5, 0.5}));
    }

    @Override
    public void run() {
        while (FEs < maxFEs) {
//            System.out.println(FEs);
            List<DoubleSolution> offspring = generateOffspring(population);
            offspring.forEach(problem::evaluate);
            population = replacement(population, offspring);
            FEs += N;
            updateEP();
        }
    }

    void updateEP() {
        for (ReferencePoint w : referencePoints) {
            if (w.members.isEmpty()) {
                continue;
            }
            List<DoubleSolution> tmp = population;
//                find the solution in the lowest level
            int lowestLevel = (int) Collections.min(w.members, Comparator.comparingInt(s -> (int) s.getAttribute("ndrank"))).getAttribute("ndrank");
            tmp = w.members.stream().filter(s -> (int) s.getAttribute("ndrank") == lowestLevel).collect(Collectors.toList());
//                if there are more than one solutions, select the closest one
            DoubleSolution bestSolution = tmp.stream().min(Comparator.comparingDouble(s -> (double) s.getAttribute("dst"))).get();
            for (int j = 0; j < EPrate.length; j++) {
                for (int i = 0; i < dim; i++) {
                    double originalX = w.evolutionPath.get(j).getVariableValue(i);
                    double newX = originalX * (1 - EPrate[j]) + bestSolution.getVariableValue(i) * EPrate[j];
                    w.evolutionPath.get(j).setVariableValue(i, newX);
                }
            }
//            add successful parameters to history
            w.members.stream().filter(s -> (boolean) s.getAttribute("isNew")).forEach(s -> {
                double[] para = (double[]) s.getAttribute("para");
                w.HistoryF[w.HistoryPtr] = para[0];
                w.HistoryCr[w.HistoryPtr] = para[1];
                w.HistoryPtr = ++w.HistoryPtr % historyLength;
            });
        }
    }

    private List<DoubleSolution> generateOffspring(List<DoubleSolution> population) {
        List<DoubleSolution> offspringPopulation = new ArrayList<>(N);
//        find the association relationship
        HashMap<DoubleSolution, ReferencePoint> AR = new HashMap<>();
        referencePoints.forEach(w -> w.members.forEach(s -> AR.put(s, w)));
        for (DoubleSolution s : population) {
//            find the corresponding reference point
            ReferencePoint w = AR.get(s);
//            determine the mating pool type
            List<DoubleSolution> matingPool = w.members;
            if (Tools.getRndGen().nextUniform(0, 1) > 0.9) matingPool = population;
//            find two neighbors
            DoubleSolution[] r12 = new DoubleSolution[2];
            r12[0] = (DoubleSolution) Tools.getRndGen().nextSample(matingPool, 1)[0];
            r12[1] = (DoubleSolution) Tools.getRndGen().nextSample(matingPool, 1)[0];
//            get parameter
            double[] para = generateParametersForDE(w);
            double F = para[0], Cr = para[1];
//            DoubleSolution future = generateFutureSolution(w, 2.0 - F); // M-0
            DoubleSolution future = generateFutureSolution(w, (1 - F) * (1 - 1.0 * FEs / maxFEs) + 1.01);
            DoubleSolution newSolution = MOEAD_EP.DE_operation(future, s, r12[0], r12[1], F, Cr);
            mutationOperator.execute(newSolution);
            newSolution.setAttribute("para", para);
            offspringPopulation.add(newSolution);
        }
        this.population.forEach(s -> s.setAttribute("isNew", false));
        offspringPopulation.forEach(s -> s.setAttribute("isNew", true));
        return offspringPopulation;
    }

    public double[] generateParametersForDE(ReferencePoint w) {
        int rndIdx = Tools.getRndGen().nextInt(0, historyLength - 1);
        double mean_f = w.HistoryF[rndIdx];
        double mean_cr = w.HistoryCr[rndIdx];
        double f, cr;
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

    private List<DoubleSolution> replacement(List<DoubleSolution> population, List<DoubleSolution> offspringPopulation) {
        List<DoubleSolution> jointPopulation = new ArrayList<>();
        jointPopulation.addAll(population);
        jointPopulation.addAll(offspringPopulation);
        Ranking<DoubleSolution> ranking = new DominanceRanking().computeRanking(jointPopulation);
        List<DoubleSolution> pop = new ArrayList<>();
        int rankingIndex = 0;
        List<DoubleSolution> candidateSet = new ArrayList<>();
        while (candidateSet.size() < N) {
            for (DoubleSolution s : ranking.getSubfront(rankingIndex)) {
                candidateSet.add(s);
                s.setAttribute("ndrank", rankingIndex);
            }
            if (candidateSet.size() <= N)
                ranking.getSubfront(rankingIndex).forEach(pop::add);
            rankingIndex++;
        }
        int lastLevelIdx = (int) candidateSet.get(candidateSet.size() - 1).getAttribute("ndrank");
        pop.addAll(selectFromLastFront(N - pop.size(), candidateSet, lastLevelIdx));
        return pop;
    }

    private void findZStarAndZNad(List<DoubleSolution> candidateSet) {
//        ideal point
        double[] z = new double[M];
        Arrays.fill(z, Double.MAX_VALUE);
        for (DoubleSolution s : candidateSet)
            for (int i = 0; i < M; i++)
                z[i] = Math.min(z[i], s.getObjective(i));
//        find extreme points
        List<DoubleSolution> extremePoints = new ArrayList<>();
        for (int i = 0; i < M; i++) {
            int finalI = i;
            DoubleSolution minSolution = Collections.min(candidateSet, Comparator.comparingDouble(a -> ASF(a, z, finalI)));
            if (extremePoints.contains(minSolution)) {
                break;
            }
            extremePoints.add(minSolution);
        }
//       find interceptions
        double[] a = null;
        if (extremePoints.size() == M) {
//        construct hyperplane
            RealMatrix E = MatrixUtils.createRealMatrix(M, M);
            for (int i = 0; i < M; i++) {
                for (int j = 0; j < M; j++) {
                    E.setEntry(i, j, extremePoints.get(i).getObjective(j));
                }
            }
            ArrayRealVector b = new ArrayRealVector(M, 1);
            DecompositionSolver solver = new LUDecomposition(E).getSolver();
            try {
                a = solver.solve(b).toArray();
                for (int i = 0; i < M; i++) {
                    if (a[i] <= 1e-10) {
                        a = null;
                        break;
                    }
                    a[i] = 1.0 / a[i];
                    if (a[i] - z[i] < 1e-10) {
                        a = null;
                        break;
                    }
                }
            } catch (Exception e) {
                a = null;
            }
        }
//        fail to construct hyperplane, use the max value instead
        if (a == null) {
            a = Arrays.copyOf(z, M);
            for (DoubleSolution s : candidateSet)
                for (int i = 0; i < M; i++)
                    a[i] = Math.max(a[i], s.getObjective(i));
        }
        z_star = z;
        z_nad = a;
    }

    private double getPerpendicularDistance(DoubleSolution s, ReferencePoint w) {
        DefaultDoubleSolution ds;
        double d1 = 0, d2 = 0;
        for (int i = 0; i < M; i++) {
            d1 += w.lambda[i] * (s.getObjective(i) - z_star[i]) / (z_nad[i] - z_star[i]);
        }
        for (int i = 0; i < M; i++) {
            d2 += Math.pow((s.getObjective(i) - z_star[i]) / (z_nad[i] - z_star[i]) - w.lambda[i] * d1, 2);
        }
        return d2;
    }

    private List<DoubleSolution> selectFromLastFront(int K, List<DoubleSolution> candidateSet, int lastLevelIdx) {
//        create a new copy
        List<ReferencePoint> availableReferencePoints = new LinkedList<>();
        referencePoints.forEach(s -> {
            availableReferencePoints.add(s);
            s.members.clear();
            s.candidates.clear();
        });
//         associate
        findZStarAndZNad(candidateSet);
        for (DoubleSolution s : candidateSet) {
            double dst = Double.MAX_VALUE;
            ReferencePoint closestReferencePoint = null;
            for (ReferencePoint w : availableReferencePoints) {
                double d2 = getPerpendicularDistance(s, w);
                if (d2 < dst) {
                    dst = d2;
                    closestReferencePoint = w;
                }
            }
            s.setAttribute("dst", dst);
            if (K == 0 || (int) s.getAttribute("ndrank") < lastLevelIdx)
                closestReferencePoint.members.add(s);
            else {
                closestReferencePoint.candidates.add(s);
            }
        }
//        niching
        List<DoubleSolution> ret = new ArrayList<>(K);
        while (ret.size() < K) {
//        find reference points with minimum count
            int minCnt = Collections.min(availableReferencePoints, Comparator.comparingInt(s -> s.members.size())).members.size();
            List<ReferencePoint> j_min_set = availableReferencePoints.stream().filter(s -> s.members.size() == minCnt).collect(Collectors.toList());
//        randomly select one
            ReferencePoint j = (ReferencePoint) Tools.getRndGen().nextSample(j_min_set, 1)[0];
//        associated solutions in the last front
            if (j.candidates.size() > 0) {
                DoubleSolution selected;
                if (j.members.isEmpty()) {
                    selected = Collections.min(j.candidates, Comparator.comparingDouble(s -> (double) s.getAttribute("dst")));
                } else {
                    selected = (DoubleSolution) Tools.getRndGen().nextSample(j.candidates, 1)[0];
                }
                j.candidates.remove(selected);
                j.members.add(selected);
                ret.add(selected);
            } else {
                availableReferencePoints.remove(j);
            }
        }
        return ret;
    }

    private double ASF(DoubleSolution s, double[] z, int index) {
        double max_ratio = Double.NEGATIVE_INFINITY;
        for (int i = 0; i < s.getNumberOfObjectives(); i++) {
            double weight = (index == i) ? 1.0 : 1e-6;
            max_ratio = Math.max(max_ratio, (s.getObjective(i) - z[i]) / weight);
        }
        return max_ratio;
    }

    private DoubleSolution generateFutureSolution(ReferencePoint w, double t) {
        DoubleSolution s = (DoubleSolution) population.get(0).copy();
        double[] coff = MOEAD_EP.getCoefficient(EPrate, t);
        for (int j = 0; j < dim; j++) {
            double V = 0;
            for (int i = 0; i < EPrate.length; i++) {
                V += coff[i] * w.evolutionPath.get(i).getVariableValue(j);
            }
            s.setVariableValue(j, V);
        }
        return s;
    }

    @Override
    public List<DoubleSolution> getResult() {
        return SolutionListUtils.getNondominatedSolutions(population);
    }

    @Override
    public String getName() {
        return "NSGA3-EP";
    }

    @Override
    public String getDescription() {
        return "Nondominated Sorting Genetic Algorithm version III";
    }

}
