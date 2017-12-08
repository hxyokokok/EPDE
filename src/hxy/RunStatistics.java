package hxy;

import org.apache.commons.math3.stat.descriptive.DescriptiveStatistics;
import org.apache.commons.math3.stat.ranking.NaturalRanking;
import org.apache.commons.math3.stat.ranking.TiesStrategy;
import org.apache.poi.hssf.usermodel.HSSFWorkbook;
import org.apache.poi.hssf.util.HSSFColor;
import org.apache.poi.ss.usermodel.*;
import org.apache.poi.ss.util.CellRangeAddress;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.*;

/**
 * Created by hxy on 2017/3/26.
 */
public class RunStatistics {
    String output_dir_path = Tools.class.getResource("..").getPath() + "/results/";
    //    List<String> use_indicator = Arrays.asList("IGD");
//    List<String> use_indicator = Arrays.asList("HV0905","IGD","Delta_p");
    List<String> use_indicator = Arrays.asList("Delta_p");

    List<String> alg_dir_list;
    List<String> benchmarks;
    int maxRuns = 20;
    //    statistics
    double[][][][] experimentData;
    double[][][] mean;
    double[][][] median;
    double[][][] stdDeviation;
    double[][][] iqr;
    double[][][] max;
    double[][][] min;
    double[][][] rankingMedian;
    // wilcoxon signed test
    int[][][] WilcoxonSignedRankResult;
    static HashMap<Integer, String> WilcoxonSignedRankResultSymbols = new HashMap<>();
    //    friedman test
    double[][] friedmantResult;

    private static HashMap<String, String> isSmallerBetter = new HashMap<>();

    static {
        isSmallerBetter.put("IGD", "IGD");
        isSmallerBetter.put("Delta_p", "Delta_p");
        WilcoxonSignedRankResultSymbols.put(1, "+");
        WilcoxonSignedRankResultSymbols.put(0, "=");
        WilcoxonSignedRankResultSymbols.put(-1, "-");
    }

    public static void main(String[] args) throws Exception {
        List<String> alg_dir_list = new ArrayList<>();

        alg_dir_list.add("NSGA3-EP");
        alg_dir_list.add("NSGA3");
//        alg_dir_list.add("MOEAD-EP");
//        alg_dir_list.add("MOEAD-SBX");
//        alg_dir_list.add("MOEAD-BLX");
//        alg_dir_list.add("MOEAD-DE");
//        alg_dir_list.add("MOEAD-HOP");

//        move the target algorithm to the first column
        for (int i = 0; i < alg_dir_list.size(); i++) {
            if (alg_dir_list.get(i).startsWith("*")) {
                String tmp = alg_dir_list.get(i).substring(1);
                alg_dir_list.set(i, alg_dir_list.get(0));
                alg_dir_list.set(0, tmp);
            }
        }
        List<String> benchmarks = new ArrayList<>();
        List<MOP> tmp = new ArrayList<MOP>();
        tmp.addAll(MOP.createManyObjectiveBenchmark("WFG"));
        tmp.addAll(MOP.createManyObjectiveBenchmark("DTLZ1-4"));
        tmp.addAll(MOP.createManyObjectiveBenchmark("DTLZ5-7"));
//        tmp.addAll(MOP.createManyObjectiveBenchmark("MINUS_DTLZ1-4"));
        for (MOP mop : tmp) {
            String pname = mop.getName();
            pname = pname.substring(0, pname.lastIndexOf("_"));
            benchmarks.add(pname);
        }
        RunStatistics rs = new RunStatistics();
//        rs.use_indicator = Arrays.asList(use_indicator);
        rs.alg_dir_list = alg_dir_list;
        rs.benchmarks = benchmarks;
        rs.loadExperimentResult();
        rs.computeDataStatistics();
        rs.computeWilcoxonSignedRankTest();
        rs.computeFriedmantTest();
        rs.write2xls();
    }

    public void write2xls() throws IOException {
        Workbook wb = new HSSFWorkbook();
        CellStyle cs1 = wb.createCellStyle();
        cs1.setFillPattern(FillPatternType.SOLID_FOREGROUND);
        CellStyle cs2 = wb.createCellStyle();
        cs2.setFillPattern(FillPatternType.SOLID_FOREGROUND);
        Row row;
        for (int indicator_idx = 0; indicator_idx < use_indicator.size(); indicator_idx++) {
            Sheet ws = wb.createSheet(use_indicator.get(indicator_idx));
//            write first row
            row = ws.createRow(0);
            row.createCell(0).setCellValue("Problem");
            row.createCell(1).setCellValue("M");
            for (int alg_idx = 0; alg_idx < alg_dir_list.size(); alg_idx++) {
                row.createCell(row.getLastCellNum()).setCellValue(alg_dir_list.get(alg_idx));
            }
//            write data
//            each row(each problem)
            for (int benchmark_idx = 0; benchmark_idx < benchmarks.size(); benchmark_idx++) {
//                find the best and the second best
//                double[] ranks = new NaturalRanking().rank(median[indicator_idx][benchmark_idx]);
                double[] ranks = rankingMedian[indicator_idx][benchmark_idx];
                row = ws.createRow(benchmark_idx + 1);
//                row.createCell(0).setCellValue(benchmarks.get(benchmark_idx));
                String[] tmps = benchmarks.get(benchmark_idx).split("_");
                String problemName = tmps[0];
                int M;
                if (problemName.equals("MINUS")) {
                    problemName = problemName + "_" + tmps[1];
                    M = Integer.parseInt(tmps[2].substring(0, tmps[2].length() - 1));
                } else {
                    M = Integer.parseInt(tmps[1].substring(0, tmps[1].length() - 1));
                }

                row.createCell(0).setCellValue(problemName);
                row.createCell(1).setCellValue(M);
                for (int alg_idx = 0; alg_idx < alg_dir_list.size(); alg_idx++) {
                    String cellvalue = String.format("%5.3e", median[indicator_idx][benchmark_idx][alg_idx]);
//                    cellvalue = '$' + cellvalue + "$ ";
//                    String cellvalue = String.format("%5.3e(%2.1e)", median[indicator_idx][benchmark_idx][alg_idx], iqr[indicator_idx][benchmark_idx][alg_idx]);
                    if (alg_idx > 0)
                        cellvalue += WilcoxonSignedRankResultSymbols.get(WilcoxonSignedRankResult[indicator_idx][benchmark_idx][alg_idx - 1]);
                    Cell c = row.createCell(row.getLastCellNum());
                    c.setCellValue(cellvalue);
                    if ((int) ranks[alg_idx] == 1) {
                        cs1.setFillForegroundColor(HSSFColor.GREY_50_PERCENT.index);
                        c.setCellStyle(cs1);
                    } else if (alg_dir_list.size() > 2 && (int) ranks[alg_idx] == 2) {
                        cs2.setFillForegroundColor(HSSFColor.GREY_25_PERCENT.index);
                        c.setCellStyle(cs2);
                    }
                }
            }

//            merge cells
            for (int row_cnt = 1; row_cnt < ws.getLastRowNum(); row_cnt += 5) {
                ws.addMergedRegion(new CellRangeAddress(row_cnt, row_cnt + 4, 0, 0));
            }
//            Wilcoxon signed rank results
            if (WilcoxonSignedRankResult != null) {
                row = ws.createRow(ws.getLastRowNum() + 1);
                String tmp = WilcoxonSignedRankResultSymbols.get(1);
                tmp += "/" + WilcoxonSignedRankResultSymbols.get(0);
                tmp += "/" + WilcoxonSignedRankResultSymbols.get(-1);
                row.createCell(0).setCellValue(tmp);
                row.createCell(1);
                row.createCell(2);
                for (int alg_idx = 1; alg_idx < alg_dir_list.size(); alg_idx++) {
                    int[] tmps = new int[]{1, 0, -1};
                    for (int candidateIdx = 0; candidateIdx < 3; candidateIdx++) {
                        int candidateResult = tmps[candidateIdx];
                        int t_sum = 0;
                        for (int problem_idx = 0; problem_idx < benchmarks.size(); problem_idx++) {
                            if (WilcoxonSignedRankResult[indicator_idx][problem_idx][alg_idx - 1] == candidateResult)
                                t_sum++;
                        }
                        tmps[candidateIdx] = t_sum;
                    }
                    row.createCell(row.getLastCellNum()).setCellValue(tmps[0] + "/" + tmps[1] + "/" + tmps[2]);
                }
            }

//            friedman ranking
            if (friedmantResult != null) {
                ws.createRow(ws.getLastRowNum() + 1);
                row = ws.createRow(ws.getLastRowNum() + 1);
                row.createCell(0).setCellValue("Friedman Ranking");
                for (int alg_idx = 0; alg_idx < alg_dir_list.size(); alg_idx++) {
                    row = ws.createRow(ws.getLastRowNum() + 1);
                    row.createCell(0).setCellValue(alg_dir_list.get(alg_idx));
                    row.createCell(1).setCellValue(friedmantResult[indicator_idx][alg_idx]);
                }
            }
        }
        wb.write(new FileOutputStream(output_dir_path + "/report.xls"));
        wb.close();
    }

    public void loadExperimentResult() throws FileNotFoundException {
        experimentData = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()][maxRuns];
        for (int indicator_idx = 0; indicator_idx < use_indicator.size(); indicator_idx++) {
            for (int alg_idx = 0; alg_idx < alg_dir_list.size(); alg_idx++) {
                String alg_name = alg_dir_list.get(alg_idx);
                System.out.println("load " + alg_name);
                File alg_dir = new File(output_dir_path, alg_name);
//                HashMap<String, ArrayList<Double>> current_algorithm_result = new HashMap<>();
                for (int benchmark_idx = 0; benchmark_idx < benchmarks.size(); benchmark_idx++) {
                    String benchmark = benchmarks.get(benchmark_idx);
//                    System.out.println(alg_name + "\t" + benchmark);
                    File test_dir = alg_dir.listFiles(it -> it.getName().startsWith(benchmark))[0];
                    Scanner sin = new Scanner(new File(test_dir, use_indicator.get(indicator_idx) + ".txt"));
                    ArrayList<Double> result = new ArrayList<>();
                    for (int runNo = 0; runNo < maxRuns; runNo++) {
                        sin.nextInt(); // run No.
                        experimentData[indicator_idx][benchmark_idx][alg_idx][runNo] = sin.nextDouble();
                    }
                    sin.close();
                }
            }
        }
    }

    private void computeFriedmantTest() {
        friedmantResult = new double[use_indicator.size()][];
        for (int indicator = 0; indicator < use_indicator.size(); indicator++) {
            if (!isSmallerBetter.containsKey(use_indicator.get(indicator))) {
                double[][] tmp = new double[benchmarks.size()][alg_dir_list.size()];
                for (int i = 0; i < benchmarks.size(); i++)
                    for (int j = 0; j < alg_dir_list.size(); j++)
                        tmp[i][j] = -1 * median[indicator][i][j];
                friedmantResult[indicator] = StateTools.friedmanRanking(tmp);
            } else {
                friedmantResult[indicator] = StateTools.friedmanRanking(median[indicator]);
            }
        }
    }

    private void computeWilcoxonSignedRankTest() {
        WilcoxonSignedRankResult = new int[use_indicator.size()][benchmarks.size()][alg_dir_list.size() - 1];
        for (int indicator = 0; indicator < use_indicator.size(); indicator++) {
            for (int problem = 0; problem < benchmarks.size(); problem++) {
                for (int algorithm = 1; algorithm < alg_dir_list.size(); algorithm++) {
                    WilcoxonSignedRankTestResult result = StateTools.wilcoxonSignedRank(experimentData[indicator][problem][0], experimentData[indicator][problem][algorithm], isSmallerBetter.containsKey(use_indicator.get(indicator)));
                    if (result.pValue < 0.05) {
                        WilcoxonSignedRankResult[indicator][problem][algorithm - 1] = result.Rplus > result.Rminus ? 1 : -1;
                    } else {
                        WilcoxonSignedRankResult[indicator][problem][algorithm - 1] = 0;
                    }
                }
            }
        }
    }

    private void computeDataStatistics() {
        mean = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        median = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        stdDeviation = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        iqr = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        min = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        max = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        rankingMedian = new double[use_indicator.size()][benchmarks.size()][alg_dir_list.size()];
        for (int indicator = 0; indicator < use_indicator.size(); indicator++) {
            for (int problem = 0; problem < benchmarks.size(); problem++) {
                for (int algorithm = 0; algorithm < alg_dir_list.size(); algorithm++) {
                    double[] tmp = experimentData[indicator][problem][algorithm];
//                    Arrays.sort(tmp);
                    DescriptiveStatistics ds = new DescriptiveStatistics();
                    for (double c : tmp) {
                        ds.addValue(c);
                    }
                    mean[indicator][problem][algorithm] = ds.getMean();
                    median[indicator][problem][algorithm] = ds.getPercentile(50);
                    stdDeviation[indicator][problem][algorithm] = ds.getStandardDeviation();
                    iqr[indicator][problem][algorithm] = ds.getPercentile(75) - ds.getPercentile(25);
                    min[indicator][problem][algorithm] = ds.getMin();
                    max[indicator][problem][algorithm] = ds.getMax();
                }
                double[] tmp = Arrays.copyOf(median[indicator][problem], alg_dir_list.size());
                if (!isSmallerBetter.containsKey(use_indicator.get(indicator)))
                    for (int i = 0; i < tmp.length; i++) tmp[i] *= -1;
                rankingMedian[indicator][problem] = new NaturalRanking(TiesStrategy.AVERAGE).rank(tmp);
            }
        }
    }

}
