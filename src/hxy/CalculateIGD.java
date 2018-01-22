package hxy;

import org.uma.jmetal.qualityindicator.QualityIndicator;
import org.uma.jmetal.qualityindicator.impl.InvertedGenerationalDistance;
import org.uma.jmetal.qualityindicator.impl.hypervolume.WFGHypervolume;
import org.uma.jmetal.util.front.Front;
import org.uma.jmetal.util.front.imp.ArrayFront;
import org.uma.jmetal.util.front.util.FrontNormalizer;
import org.uma.jmetal.util.front.util.FrontUtils;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.lang.reflect.Field;
import java.util.*;

/**
 * Created by hxy on 2017/3/26.
 */
public class CalculateIGD {
    public static void main(String[] args) throws IOException, NoSuchFieldException, IllegalAccessException {
        int max_runs = 20;
        String algName = "MOEAD_EP";
        String experiment_base_dir = new File(Tools.class.getResource("..").getPath()).getParentFile().getParentFile().getParentFile().getPath()+ "/results/" + algName;

//        String experiment_base_dir = Tools.class.getResource("..").getPath() + "/results/" + algName;
        List<String> use_indicator = new ArrayList<>();
        use_indicator.add("IGD");
        for (File one_run_dir : new File(experiment_base_dir).listFiles()) {
            String[] tmp = one_run_dir.getName().split("_");
            String problem_name = tmp[0];
            int M = -1;
            if (problem_name.equals("MINUS")) {
                problem_name += '_' + tmp[1];
                M = Integer.parseInt(tmp[2].substring(0, tmp[2].length() - 1));
            } else {
                M = Integer.parseInt(tmp[1].substring(0, tmp[1].length() - 1));
            }
            HashMap<String, List<Double>> result_map = new HashMap<>();
            for (String indicator_name : use_indicator) {
                QualityIndicator quality_indicator = null;
                String pf_fname = problem_name + "_" + M + ".txt";
//                pf_fname = "/pf/" + pf_fname;
                pf_fname = new File(Tools.class.getResource("..").getPath()).getParentFile().getParentFile().getParentFile().getPath()+"/pf/" + pf_fname;
                switch (indicator_name) {
                    case "IGD":
                        quality_indicator = new InvertedGenerationalDistance(pf_fname, 1.0);
                        break;
                    case "HV":
                        System.err.println("undefined HV in java");
                        break;
                }
                result_map.put(indicator_name, new ArrayList<>());
                double[] result_list = new double[max_runs + 1];
                Arrays.fill(result_list, -1);
//                read existed indicator file
                File existed_indicator_file = new File(one_run_dir, indicator_name + ".txt");
                if (existed_indicator_file.exists()) {
                    Scanner sin = new Scanner(existed_indicator_file);
                    while (sin.hasNext()) {
                        int runNo = sin.nextInt();
                        result_list[runNo] = sin.nextDouble();
                    }
                    sin.close();
                }
                for (int runs = 1; runs <= max_runs; runs++) {
                    System.out.print(problem_name);
                    System.out.print("\tM" + M);
                    System.out.print("\t@" + indicator_name + "\t#" + runs);
                    double qi;
                    if (result_list[runs] != -1) {
                        System.out.print("\t[skipped]");
                        qi = result_list[runs];
                    } else {
                        File currentFile = new File(one_run_dir, "FUN" + runs + ".tsv");
                        List solutionList = FrontUtils.convertFrontToSolutionList(new ArrayFront(currentFile.getPath()));
                        qi = (double) quality_indicator.evaluate(solutionList);
                    }
                    System.out.println("\t:" + qi);
                    result_list[runs] = qi;
                }
//                write result to files
                if (!existed_indicator_file.exists())
                    existed_indicator_file.createNewFile();
                PrintWriter out = new PrintWriter(existed_indicator_file);
                for (int runs = 1; runs <= max_runs; runs++) {
                    out.println(runs + "\t" + result_list[runs]);
                }
                out.close();
                System.out.println("writing result to [" + existed_indicator_file.getPath() + "]");
            }
        }
    }
}
