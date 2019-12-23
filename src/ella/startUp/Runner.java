package ella.startUp;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.indexing.IndexCreator;
import ella.model.io.MyParameters;

import java.io.File;

public class Runner {

    public static void run(String[] args) {
        String mode = args[0];
        switch (mode) {
            case "makedb":
                if (!IndexOptionHandler.run(args))
                    System.exit(1);
                File inFile = IndexOptionHandler.inFile;
                File edbFile = IndexOptionHandler.dbFile;
                int indexCores = IndexOptionHandler.cores;
                Long size = IndexOptionHandler.size;
                int[] min = IndexOptionHandler.min;
                int p = min != null ? min[0] : 1;
                int q = min != null ? min[1] : 1;
                new IndexCreator().run(inFile, edbFile, 3, MyParameters.NO_SEED_SHAPE, p, q, indexCores, size);
                break;
            case "blastx":
                if (!AlignOptionHandler.run(args, true))
                    System.exit(1);
                File indexFile = AlignOptionHandler.indexFile;
                File queryFile = AlignOptionHandler.queryFile;
                MyParameters.FRAMESHIFT_PENALTY = MyParameters.FRAMESHIFT_PENALTY == 0 ? null : MyParameters.FRAMESHIFT_PENALTY;
                int m = AlignOptionHandler.m;
                int alignCores = AlignOptionHandler.cores;
                File daaFile = AlignOptionHandler.outputFile;
                new Aligner(indexFile, queryFile, daaFile, m, alignCores).run();
                break;
            default:
                System.err.println("Wrong mode: " + mode);
        }
    }

}
