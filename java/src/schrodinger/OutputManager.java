package schrodinger;

import java.io.IOException;
import java.io.PrintWriter;
import java.io.Writer;
import java.nio.charset.StandardCharsets;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.nio.file.StandardOpenOption;
import java.time.LocalDateTime;
import java.time.format.DateTimeFormatter;
import java.util.Set;
import java.util.concurrent.ConcurrentHashMap;

/**
 * Orchestrates writing to a common output file and per‑level files.
 * Uses buffered writing for performance; relies on JVM shutdown for final flush.
 */
public class OutputManager {
    // Common output (static)
    private static PrintWriter commonWriter;
    private static final DateTimeFormatter TIMESTAMP = DateTimeFormatter.ofPattern("yyyy-MM-dd HH:mm:ss.SSS");

    // Registry of all per‑level writers for shutdown cleanup
    private static final Set<LevelWriter> levelWriters = ConcurrentHashMap.newKeySet();

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
     * @param filePath destination file (overwritten); if null, common output is disabled.
     */
    public static synchronized void initCommonOutputFile(String filePath) {
        closeCommon();
        if (filePath != null) {
            try {
                Path path = Paths.get(filePath);
                // Ensure parent directories exist if they don't already
                if (path.getParent() != null) {
                    Files.createDirectories(path.getParent());
                }

                Writer writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE, StandardOpenOption.TRUNCATE_EXISTING, StandardOpenOption.WRITE);
                commonWriter = new PrintWriter(writer);
            } catch (IOException e) {
                throw new RuntimeException(e);
            }
        }
    }

    private static synchronized void closeCommon() {
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
        }
    }

    public static synchronized void writeBlankLine() {
        write("");
    }

    /**
     * Write a raw line (no timestamp) to the common output – useful for CSV data.
     * Thread‑safe.
     */
    public static synchronized void writeData(String line) {
        if (commonWriter != null) {
            commonWriter.println(line);
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
                Path path = Paths.get(filePath);
                if (path.getParent() != null) {
                    Files.createDirectories(path.getParent());
                }

                Writer writer = Files.newBufferedWriter(path, StandardCharsets.UTF_8,
                        StandardOpenOption.CREATE, StandardOpenOption.APPEND, StandardOpenOption.WRITE);
                this.writer = new PrintWriter(writer);
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
        }

        /**
         * Write a raw line (no timestamp) to this level's file.
         */
        public void raw(String line) {
            writer.println(line);
        }

        /**
         * Closes the writer and removes it from the shutdown registry.
         * Call this if you are done with a level writer before the JVM exits.
         */
        public void close() {
            writer.flush();
            writer.close();
            levelWriters.remove(this);
        }
    }
}