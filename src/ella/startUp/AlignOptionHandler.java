package ella.startUp;

import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.io.MyParameters;

import java.io.File;
import java.util.ArrayList;

public class AlignOptionHandler {

    public static File indexFile, queryFile, outputFile;
    public static int cores = Runtime.getRuntime().availableProcessors();
    public static int m = 10;
    public static String matrix = "BLOSUM80";
    private static Integer gop = null, gep = null;

    private static Integer F = null, minCoverage = null;
    private static Double eValue = null;

    public static boolean run(String[] args, boolean setParameters) {

        if (args.length == 1 && (args[1] == "help" || args[1] == "-h" || args[1] == "-help")) {
            printOptions();
            return false;
        }

        for (int i = 1; i < args.length; i++) {
            String option = args[i];
            switch (option) {
                case "--threads":
                case "-p":
                    try {
                        cores = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -p: not an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--out":
                case "-o":
                    try {
                        outputFile = new File(args[i + 1]);
                        outputFile.createNewFile();
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -o: could not create output file " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--db":
                case "-d":
                    try {
                        indexFile = new File(args[i + 1]);
                        if (!indexFile.exists()) {
                            System.err.print("ERROR parsing option -d: file does not exist " + (args[i + 1]));
                            return false;
                        }
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -d: not a file " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--in":
                case "-i":
                    try {
                        queryFile = new File(args[i + 1]);
                        if (!queryFile.exists()) {
                            System.err.print("ERROR parsing option -i: file does not exist " + (args[i + 1]));
                            return false;
                        }
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -i: not a file " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "-m":
                    try {
                        m = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -m: not an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "-F":
                    try {
                        F = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -F: not an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--evalue":
                case "-e":
                    try {
                        eValue = Double.parseDouble(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -e: not a double " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--coverage":
                case "-c":
                    try {
                        minCoverage = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option -c: not an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--matrix":
                    try {
                        matrix = args[i + 1];
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option --matrix: no type specified " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--gop":
                    try {
                        gop = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option --gop: no an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                case "--gep":
                    try {
                        gep = Integer.parseInt(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR parsing option --gep: no an integer " + (args[i + 1]));
                    }
                    i++;
                    break;
                default:
                    System.out.println("ERROR: Unknown option " + option);
                    System.exit(1);
            }
        }

        if (!BlastStatisticsHelper.matrix2penalties.containsKey(matrix)) {
            System.err.println("ERROR: unknown scoring matrix (Option --matrix)");
            return printMatrixSupport();
        }

        if (gop == null)
            gop = getDefaultGOP();
        if (gep == null)
            gep = getDefaultGEP();
        if (!checkPenalties()) {
            System.err.println("ERROR: unsupported gap-open pr gap-extend penalties (Option --gop/--gep)");
            return printMatrixSupport();
        }
        if (indexFile == null) {
            System.err.println("ERROR: missing path DIAMOND database file (Option --db/-d)");
            return printOptions();
        }
        if (queryFile == null) {
            System.err.println("ERROR: missing query input file in FASTA or FASTQ format (Option --in/-i)");
            return printOptions();
        }
        if (queryFile == null) {
            System.err.println("ERROR: missing path to output file (Option --out/-o)");
            return printOptions();
        }

        MyParameters.CPUS = cores;
        if (setParameters) {
            MyParameters.setScoringMatrix(matrix, gop, gep);
            if (F != null)
                MyParameters.FRAMESHIFT_PENALTY = F;
            if (eValue != null)
                MyParameters.MAX_EVALUE = eValue;
            if (minCoverage != null)
                MyParameters.MIN_COVERAGE = minCoverage;
        }

        return true;

    }

    private static boolean checkPenalties() {
        ArrayList<String> penalties = BlastStatisticsHelper.matrix2penalties.get(matrix);
        for (String p : penalties) {
            if (p.equals(gop + "/" + gep))
                return true;
        }
        return false;
    }

    private static Integer getDefaultGEP() {
        String s = BlastStatisticsHelper.matrix2penalties.get(matrix).get(0);
        s = s.split("/")[1];
        return Integer.parseInt(s);
    }

    private static Integer getDefaultGOP() {
        String s = BlastStatisticsHelper.matrix2penalties.get(matrix).get(0);
        s = s.split("/")[0];
        return Integer.parseInt(s);
    }

    private static boolean printOptions() {
        System.out.println("Mandatory Options:");
        System.out.println("--db/-d \t path to the indexed protein reference database file in EDB format");
        System.out.println("--in/-i \t path to the query read file in FASTA/FASTQ format");
        System.out.println("--out/-o \t path of the reported output file in DAA format");
        System.out.println("Optional");
        System.out.println("--threads/-p \t number of parallel threads (default: #cpus)");
        System.out.println("-F \t penalty for switching frames during blastx alignments (default: 15)");
        System.out.println("-m \t maximum initial matches per query position (default: 10)");
        System.out.println("-e \t maximum e-value for reported alignments (default: 0.001)");
        System.out.println("-c \t minimum percentage of reference positions covered by reported alignments (default: 0)");
        System.out.println("--matrix \t scoring matrix (default: BLOSUM80)");
        System.out.println("--gop \t gap open penalty");
        System.out.println("--gep \t gap extension penalty");
        return false;
    }

    private static boolean printMatrixSupport() {
        System.out.println("The following scoring matrices are supported:");
        System.out.println("Matrix\tPenalties(gap-open/gap-extend)");
        System.out.println("---");
        for (String key : BlastStatisticsHelper.matrix2penalties.keySet()) {
            StringBuilder s = new StringBuilder(key + "\t");
            for (String p : BlastStatisticsHelper.matrix2penalties.get(key))
                s.append(p + " ");
            System.out.println(s);
        }
        return false;
    }

}
