package hxy;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.stat.inference.WilcoxonSignedRankTest;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;

/**
 * Created by hxy on 2017/4/10.
 */
public class StateTools {
    private static NaturalRanking ranking = new NaturalRanking(TiesStrategy.AVERAGE);

    public static double[] friedmanRanking(double[][] data) {
        int n = data.length; // problem number
        int k = data[0].length; // algorithm number

        double[] meanRanking = new double[k];
        double[][] r = new double[n][];
        for (int problem = 0; problem < n; problem++) {
            r[problem] = ranking.rank(data[problem]);
        }
        for (int algorithm = 0; algorithm < k; algorithm++) {
            for (int problem = 0; problem < n; problem++) {
                meanRanking[algorithm] += r[problem][algorithm];
            }
            meanRanking[algorithm] /= n;
        }
        return meanRanking;
    }

    public static WilcoxonSignedRankTestResult wilcoxonSignedRank(double[] x, double[] y, boolean isSmallerBetter) {
        WilcoxonSignedRankTestResult result = new WilcoxonSignedRankTestResult();
        if (x.length != y.length) {
            try {
                throw new Exception("x.length ~= y.length");
            } catch (Exception e) {
                e.printStackTrace();
            }
        }
        double[] diff = new double[x.length];
        for (int i = 0; i < diff.length; i++)
            diff[i] = x[i] - y[i];
        double[] diffAbs = new double[x.length];
        for (int i = 0; i < diff.length; i++)
            diffAbs[i] = Math.abs(diff[i]);
        double[] ranks = ranking.rank(diffAbs);
        double Wplus = 0.0D;
        for (int i = 0; i < diff.length; ++i) {
            if (diff[i] > 0.0D) {
                Wplus += ranks[i];
            } else if (diff[i] == 0) {
                Wplus += ranks[i] / 2;
            }
        }
        int N = x.length;
        double Wminus = (double) (N * (N + 1)) / 2.0D - Wplus;
        double T = Math.min(Wminus, Wplus);
        double ES = (double) (N * (N + 1)) / 4.0D;
        double VarS = ES * ((double) (2 * N + 1) / 6.0D);
        double z_value = (T - ES - 0.5D) / Math.sqrt(VarS);
        NormalDistribution standardNormal = new NormalDistribution(0.0D, 1.0D);
        if (isSmallerBetter) {
            result.Rplus = Wminus;
            result.Rminus = Wplus;
        } else {
            result.Rplus = Wplus;
            result.Rminus = Wminus;
        }
        result.pValue = 2.0D * standardNormal.cumulativeProbability(z_value);
        return result;
    }
}

class WilcoxonSignedRankTestResult {
    double Rplus, Rminus;
    double pValue;
}
