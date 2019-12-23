package ella.model.aligner.aligning.seeding;

import java.util.ArrayList;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.aligning.query.QueryContainer;
import ella.model.aligner.indexing.IndexText;
import ella.model.aligner.indexing.SuffixArrayTable;
import ella.model.aligner.utils.CodonTranslator;
import ella.model.aligner.utils.MyFrame;
import ella.model.aligner.utils.MyStrand;
import ella.model.aligner.utils.suffixArray.ImprovedBinarySearch;
import ella.model.aligner.utils.wrapper.ObjectArrayWrapper;

public class BasicSeedFinder {

    private static String query = "Escherichia-coli-str.-K-12-substr.-MG1655,-synthetic-genome_4174425_aligned_840_R_0_12711_0";
    private static String reference = "NP_418414.1";
    private boolean debug = false;

    private final MatchExtender matchExtender = new MatchExtender();
    private final ImprovedBinarySearch search = new ImprovedBinarySearch();

    public void run(QueryContainer qC, IndexText indexText, SuffixArrayTable saTable, int m, Aligner aligner) {

        int sepDepth = saTable.getSeparationDepth();
        int len = 100; // saTable.getSuffixLength();
        String dnaPlus = qC.getStrandSpecificDNA(MyStrand.FWD);
        String dnaMinus = qC.getStrandSpecificDNA(MyStrand.REV);

        for (MyFrame frame : MyFrame.values()) {

            String dna = frame.getNumericalID() > 0 ? dnaPlus : dnaMinus;
            String aa = CodonTranslator.translate(dna, frame.getNumericalID());
            DiagonalTable diagonalTable = new DiagonalTable();

            for (int i = 0; i < aa.length() - sepDepth; i += 1) {

                // finding all matching positions for k-mer at position i
                int j = i + len < aa.length() ? i + len : aa.length();
                String kMer = aa.substring(i, j);

                ArrayList<Object[]> matches = saTable.findPattern(kMer, 1, m, false, qC, aligner, search);
                aligner.counter += matches.size();

                // storing all matches for position i
                int maxMatchLength = 0;
                for (Object[] match : matches) {

                    // checking if coordinate has already been covered by previous seed extensions
                    if (diagonalTable.coversCoordinate(i, (long) match[0]))
                        continue;
                    Object[] extResult = matchExtender.run(match, indexText, i, aa, frame, m, aligner, diagonalTable);

                    if (extResult == null)
                        continue;
                    Seed seed = (Seed) extResult[1];

                    diagonalTable.addCoordinate(seed.getQueryEnd() / 3, seed.getIndexEnd());
                    ObjectArrayWrapper info = (ObjectArrayWrapper) extResult[0];
                    qC.addSeed(seed, info, frame.getStrand());
                    maxMatchLength = maxMatchLength < (int) match[1] ? (int) match[1] : maxMatchLength;

                }

                i += maxMatchLength;

            }

        }

    }

}
