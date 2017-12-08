package hxy;

import org.uma.jmetal.algorithm.Algorithm;
import org.uma.jmetal.solution.DoubleSolution;
import org.uma.jmetal.util.fileoutput.SolutionListOutput;
import org.uma.jmetal.util.fileoutput.impl.DefaultFileOutputContext;

import java.io.File;
import java.text.SimpleDateFormat;
import java.util.*;
import java.util.stream.IntStream;

/**
 * Created by hxy on 2017/3/27.
 */
public class RunExperiment {
    public static void save(Algorithm<List<DoubleSolution>> alg, String experiment_dir_path, String alg_name, MOP p, int runNo) {
        List<DoubleSolution> result = alg.getResult();
        File test_dir = new File(experiment_dir_path + "/" + alg_name + "/" + p.getName());
        if (!test_dir.exists()) {
            test_dir.mkdirs();
        }
        new SolutionListOutput(result).setSeparator("\t").setFunFileOutputContext(new DefaultFileOutputContext(test_dir.getPath() + "/FUN" + runNo + ".tsv")).print();
        System.out.printf("%s on %s #%d at %s\n", alg_name, p.getName(), runNo, new SimpleDateFormat("yyyy-MM-dd HH:mm:ss").format(new Date()));
    }

    public static void main(String[] args) throws Exception {
        int maxRuns = 20;
        List<MOP> plist = new ArrayList<>();
//        plist.addAll(MOP.createManyObjectiveBenchmark("WFG"));
        plist.addAll(MOP.createManyObjectiveBenchmark("DTLZ1-4"));
//        plist.addAll(MOP.createManyObjectiveBenchmark("DTLZ5-7"));
//        plist.addAll(MOP.createManyObjectiveBenchmark("MINUS_DTLZ1-4"));

        String experiment_dir_path = Tools.class.getResource("..").getPath() + "/results/";

        IntStream.range(1, maxRuns + 1).parallel().forEach((int runNo) -> {
            for (MOP p : plist) {
                String alg_name = "MOEAD_EP";
                Algorithm alg = null;
                boolean isSuccess;
                do {
                    try {
                        isSuccess = true;
                        alg = new MOEAD_EP(p);
                        alg.run();
                    } catch (Exception e) {
                        isSuccess = false;
                        System.out.println("!!!error!!!");
                    }
                } while (!isSuccess);
                save(alg, experiment_dir_path, alg_name, p, runNo);
            }
        });

//        Runtime.getRuntime().exec("shutdown -s -f -t 120");
    }
}
