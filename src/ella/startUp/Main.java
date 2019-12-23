package ella.startUp;

public class Main {

    public static void main(String[] args) {
        if (args.length == 0)
            printModes("ERROR: No mode specified. The available modes are:");
        String mode = args[0];
        if (!mode.equals("makedb") && !mode.equals("blastx"))
            printModes("ERROR: Unknown mode " + args[0] + ". The available modes are:");
        else
            Runner.run(args);
        System.exit(0);
    }

    public static void printModes(String message) {
        System.out.println(message);
        System.out.println("makedb\t use this mode in order to create an index for a protein reference database");
        System.out.println("blastx\t use this mode to align a dna query FASTA file against an already constructed index");
        System.exit(1);
    }

}
