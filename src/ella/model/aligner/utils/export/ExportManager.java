package ella.model.aligner.utils.export;

import ella.model.aligner.utils.daaMerger.DaaHit;
import ella.model.aligner.utils.daaMerger.ReadHits;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.ConcurrentHashMap;

public class ExportManager {

    public enum FORMAT {BLAST_TAB}

    public static void run(File out, File daaFile, FORMAT format) {
        new Thread(() -> {
            switch (format) {
                case BLAST_TAB:
                    runBlastTab(out, daaFile);
                    break;
            }
        }).start();
    }

    private static void runBlastTab(File daaFile, File out) {
//        DAA_Reader daaReader = new DAA_Reader(daaFile, false);
//        DAA_Header daaHeader = daaReader.getDAAHeader();
//        AlignmentStatisticsHelper.init(MyParameters.SCORING_MATRIX.getType(), daaHeader.getGapOpen(), daaHeader.getGapExtend(), daaHeader.getDbSeqsUsed());
//        try {
//            RandomAccessFile raf = new RandomAccessFile(daaFile, "r");
//            BlastTab.run(out, extractReads(daaReader.parseAllHits(AlignOptionHandler.cores));
//            raf.close();
//        } catch (Exception e) {
//            e.printStackTrace();
//        }
    }

    private static ArrayList<DaaHit> extractReads(ConcurrentHashMap<String, ReadHits> readHits) {
//        ArrayList<DAA_Hit> allHits = new ArrayList<>();
//        for (String read : readHits.keySet()) {
//            readHits.get(read).getAllHits();
//        }
        return null;
    }


}
