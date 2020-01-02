package ella.model.io;

import java.io.*;
import java.util.zip.GZIPInputStream;
import java.util.zip.ZipException;

public class FastaIterator {

    private byte[] buffer = new byte[1024];
    private int breakPoint = -1, readChars = 0;
    private InputStream is;

    private StringBuffer accBuffer = new StringBuffer(), sequenceBuffer = new StringBuffer();

    public FastaIterator(File file) {
        try {
            try {
                is = new BufferedInputStream(new GZIPInputStream(new FileInputStream(file)));
            } catch (ZipException e) {
                is = new BufferedInputStream(new FileInputStream(file));
            }
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public StringBuffer[] next() {
        try {

            accBuffer = new StringBuffer();
            sequenceBuffer = new StringBuffer();

            // scanning entire file
            boolean readingSequence = true, indexHeader = false;
            while (breakPoint >= 0 || (readChars = is.read(buffer)) != -1) {
                for (int i = Math.max(breakPoint, 0); i < readChars; ++i) {
                    char c = (char) buffer[i];
                    if (c == '>') {
                        readingSequence = false;
                        indexHeader = true;
                        if (accBuffer.length() > 0) {
                            StringBuffer[] entry = {accBuffer, sequenceBuffer};
                            breakPoint = i;
                            return entry;
                        }
                    } else if (readingSequence) {
                        if (c != '\n')
                            sequenceBuffer.append(c);
                    } else if (c == '\n') {
                        readingSequence = true;
                        indexHeader = false;
                    } else if (c == ' ')
                        indexHeader = false;
                    else if (indexHeader) {
                        accBuffer.append(c);
                    }
                }

                breakPoint = -1;

            }

            if (readChars == -1 && accBuffer.length() > 0 && sequenceBuffer.length() > 0) {
                StringBuffer[] entry = {accBuffer, sequenceBuffer};
                return entry;
            }

        } catch (Exception e) {
            e.printStackTrace();
        }

        return null;

    }

    public boolean hasNext() {
        return readChars != -1;
    }

    public void close() {
        if (is != null) {
            try {
                is.close();
            } catch (IOException e) {
                e.printStackTrace();
            }
        }
    }

}
