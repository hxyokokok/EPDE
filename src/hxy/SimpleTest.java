package hxy;

import org.uma.jmetal.qualityindicator.impl.InvertedGenerationalDistance;
import org.uma.jmetal.solution.DoubleSolution;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import static java.lang.System.*;

/**
 * Created by hxy on 2017/3/27.
 */
public class SimpleTest {
    public static void main(String[] args) throws Exception {
        List<MOP> plist = new ArrayList<>();
        plist.addAll(MOP.createManyObjectiveBenchmark("WFG"));
//        plist.addAll(MOP.createManyObjectiveBenchmark("DTLZ1-4"));
//        plist.addAll(MOP.createManyObjectiveBenchmark("DTLZ5-7"));
//        plist.addAll(MOP.createManyObjectiveBenchmark("MINUS_DTLZ1-4"));
        for (MOP p : plist) {
//            NSGA3_EP alg = new NSGA3_EP(p);
            MOEAD_EP alg = new MOEAD_EP(p);
            long start = currentTimeMillis();
            alg.run();
            double elapsedTime = (currentTimeMillis() - start) / 1000.0;
            List<DoubleSolution> result = alg.getResult();
            String pf_name = new File(Tools.class.getResource("..").getPath()).getParentFile().getParentFile().getParentFile().getPath()+"/pf/" + p.getPfName();
//            String pf_name = "/pf/" + p.getPfName();
//            config IGD
            InvertedGenerationalDistance<DoubleSolution> IGD_measure = new InvertedGenerationalDistance<DoubleSolution>(pf_name, 1);
            double IGD = IGD_measure.evaluate(result);
            out.println(alg.getName() + "\t" + p.getName() + "\tIGD:" + IGD + "\tTime:" + elapsedTime + "s");
        }
    }
}
