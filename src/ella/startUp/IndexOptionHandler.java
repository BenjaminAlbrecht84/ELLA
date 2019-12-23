package ella.startUp;

import java.io.File;

public class IndexOptionHandler {

    public static File inFile, dbFile;
    public static int cores = Runtime.getRuntime().availableProcessors();
    public static Long size = null;
    public static int[] min = null;

    public static boolean run(String[] args) {

        if (args.length == 2 && (args[1] == "help" || args[1] == "-h" || args[1] == "-help")) {
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
                        System.err.print("ERROR: not an integer " + (args[i + 1])+" (Option --threads/-p)");
                    }
                    i++;
                    break;
                case "--db":
                case "-d":
                    try {
                        dbFile = new File(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR: not a file " + (args[i + 1]) + " (Option --db/-d)");
                    }
                    i++;
                    break;
                case "--in":
                case "-i":
                    try {
                        inFile = new File(args[i + 1]);
                    } catch (Exception e) {
                        System.err.print("ERROR: not a file " + (args[i + 1]) + " (Option --in/-i)");
                    }
                    i++;
                    break;
                case "--size":
                case "-s":
                    size = parseSize(args[i + 1]);
                    i++;
                    break;
                case "--minimizer":
                case "-m":
                    min = parseMinimizer(args[i + 1]);
                    i++;
                    break;
                default:
                    System.out.println("ERROR: Unknown option " + option);
                    return false;
            }
        }

        if (dbFile == null) {
            System.err.println("ERROR: missing path to indexed database file (Option --db/-d)");
            printOptions();
            return false;
        }
        if (inFile == null) {
            System.err.println("ERROR: missing query input file in FASTA or FASTQ format (Option --query/-q)");
            printOptions();
            return false;
        }

        return true;
    }

    private static int[] parseMinimizer(String s) {
        String[] pq = s.substring(1, s.length() - 1).split(",");
        try {
            int p = Integer.parseInt(pq[0]);
            int q = Integer.parseInt(pq[1]);
            int[] result = {p, q};
            return result;
        } catch (Exception e) {
            System.err.print("ERROR: not an integer " + pq[0] + " " + pq[1]);
            printOptions();
        }
        return null;
    }

    private static Long parseSize(String s) {
        long factor = 1;
        char c = s.charAt(s.length() - 1);
        switch (c) {
            case 'g':
                factor *= Math.pow(10, 9);
                break;
            case 'm':
                factor *= Math.pow(10, 6);
                break;
            default:
                System.err.println("ERROR: unknown volume size delimiter " + c + " (Option --size/-s)");
                printOptions();
        }

        try {
            long size = Integer.parseInt(s.substring(0, s.length() - 1));
            size *= factor;
            return size;
        } catch (Exception e) {
            System.err.print("ERROR: not an integer " + (s.substring(0, s.length() - 1)) + " (Option --size/-s)");
            printOptions();
        }
        return null;
    }

    private static void printOptions() {
        System.out.println("Mandatory Options:");
        System.out.println("--in/-i \t path to the input protein reference database file in FASTA format");
        System.out.println("--db/-d \t path to the output indexed reference database file in EDB format");
        System.out.println("--size/-s \t maximal size of each volume; use delimiter 'g' (gb) or 'm' (mb)");
        System.out.println("--threads/-p \t number of parallel threads (default: #cpus)");
    }

}
