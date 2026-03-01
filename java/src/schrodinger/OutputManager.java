package schrodinger;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Orchestrates writing to a common output file and per‑level files.
 * All files are flushed after each write (adjust if performance becomes an issue).
 */
public class OutputManager {
    // Common output (static)
    private static PrintWriter commonWriter;
    private static final DateTimeFormatter TIMESTAMP =
            DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss");

    // Registry of all per‑level writers for shutdown cleanup
    private static final Set<LevelWriter> levelWriters =
            ConcurrentHashMap.newKeySet();

    static {
        Runtime.getRuntime().addShutdownHook(new Thread(() -> {
            closeCommon();
            for (LevelWriter lw : levelWriters) {
                lw.close();
            }
        }));
    }

    private OutputManager() {  // prevent instantiation
    }

    /**
     * Initializes the common output file.
     *
     * @param filePath destination file (appended); if null, common output is disabled.
     */
    public static void initCommon(String filePath) {
        closeCommon();
        if (filePath != null) {
            try {
                commonWriter = new PrintWriter(new FileWriter(filePath, false));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

    }

    private static void closeCommon() {
        if (commonWriter != null) {
            commonWriter.flush();
            commonWriter.close();
            commonWriter = null;
        }
    }

    /**
     * Write a timestamped message to the common output.
     * Thread‑safe (synchronized).
     */
    public static synchronized void write(String msg) {
        if (commonWriter != null) {
            commonWriter.println(timestamp() + " " + msg);
            commonWriter.flush();
        }
    }

    /**
     * Write a raw line (no timestamp) to the common output – useful for CSV data.
     * Thread‑safe.
     */
    public static synchronized void writeData(String line) {
        if (commonWriter != null) {
            commonWriter.println(line);
            commonWriter.flush();
        }
    }

    // --- Per‑level writer factory ---

    /**
     * Creates a new writer for a specific energy level.
     *
     * @param filePath destination file (appended)
     */
    public static LevelWriter createLevelWriter(String filePath) {
        LevelWriter lw = new LevelWriter(filePath);
        levelWriters.add(lw);
        return lw;
    }

    // --- Helper ---
    private static String timestamp() {
        return LocalDateTime.now().format(TIMESTAMP);
    }

    // --- Per‑level writer class ---
    public static class LevelWriter {
        private final PrintWriter writer;

        private LevelWriter(String filePath) {
            try {
                this.writer = new PrintWriter(new FileWriter(filePath, true));
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }

        /**
         * Write a timestamped message to this level's file.
         * Not synchronized – assume single‑threaded use per level.
         */
        public void log(String msg) {
            writer.println(timestamp() + " " + msg);
            writer.flush();
        }

        /**
         * Write a raw line (no timestamp) to this level's file.
         */
        public void raw(String line) {
            writer.println(line);
            writer.flush();
        }

        // Called by shutdown hook
        private void close() {
            writer.flush();
            writer.close();
        }
    }
}