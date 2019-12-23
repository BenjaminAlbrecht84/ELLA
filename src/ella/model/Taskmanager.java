package ella.model;

import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.indexing.IndexCreator;
import ella.model.io.MyParameters;
import ella.model.utils.PausableThreadPoolExecutor;
import ella.model.utils.ProfileThread;
import ella.model.utils.RuntimeThread;
import ella.presenter.Presenter;
import ella.startUp.AlignOptionHandler;
import ella.startUp.IndexOptionHandler;
import javafx.application.Platform;
import javafx.beans.property.*;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;

import java.io.File;
import java.util.ArrayList;
import java.util.concurrent.Executors;
import java.util.concurrent.ThreadPoolExecutor;

public class Taskmanager {

    private PausableThreadPoolExecutor executor = new PausableThreadPoolExecutor((ThreadPoolExecutor) Executors.newFixedThreadPool(1));

    private StringProperty status = new SimpleStringProperty();
    private int idCounter = 1;
    private DoubleProperty runningTaskProgress = new SimpleDoubleProperty(0.);
    private EllaTask runningTask = null;
    private ObservableList<EllaTask> taskQueue = FXCollections.observableArrayList();
    private BooleanProperty isRunning = new SimpleBooleanProperty(false);

    private Presenter presenter;

    public void pauseResume() {
        if (executor.isPaused()) {
            executor.resume();
            for (EllaTask task : taskQueue)
                task.start();
            isRunning.set(true);
        } else {
            setStatus("Paused!");
            executor.pause();
            isRunning.set(false);
        }
    }

    public void stop() {
        for (EllaTask task : taskQueue)
            task.stop();
        executor.shutdown();
    }

    public void addTask(String settings, String description, String mode) {
        EllaTask task = new EllaTask(idCounter++, mode, description, settings);
        switch (mode) {
            case "blastx":
                if (AlignOptionHandler.run((mode + " " + settings).split("\\s+"), false)) {
                    taskQueue.add(task);
                    presenter.reportNewTask(task.getId(), "blastx", settings, description);
                    if (isRunning.get())
                        task.start();
                }
                break;
            case "index":
                if (IndexOptionHandler.run((mode + " " + settings).split("\\s+"))) {
                    taskQueue.add(task);
                    presenter.reportNewTask(task.getId(), "index", settings, description);
                    if (isRunning.get())
                        task.start();
                }
                break;
        }
    }

    public void cancelTasks(ArrayList<Integer> taskIDs) {
        for (int id : taskIDs) {
            EllaTask task = getTask(id);
            if (task != null) {
                task.cancel();
            }
        }
    }

    public void removeTasks(ArrayList<Integer> taskIDs) {
        for (int id : taskIDs) {
            EllaTask task = getTask(id);
            if (task != null) {
                task.stop();
                taskQueue.remove(task);
            }
        }
    }

    public EllaTask getTask(int id) {
        for (EllaTask t : taskQueue) {
            if (t.getId() == id)
                return t;
        }
        return null;
    }

    public void setStatus(String s) {
        Platform.runLater(() -> {
            if (runningTask != null) {
                presenter.reportLogInfo(runningTask.getId(), s);
                statusProperty().set("Task " + runningTask.getId() + ": " + s);
            }
        });
    }

    public void reportFinish() {
        runningTask.stop();
    }

    public int getNextTaskId() {
        return idCounter;
    }

    public void moveTasksUp(ArrayList<Integer> taskIDs) {
        EllaTask t1 = null, t2 = null;
        for (int id : taskIDs) {
            for (EllaTask t : taskQueue) {
                if (t.getId() == id) {
                    t1 = t;
                    break;
                } else
                    t2 = t;
            }
            if (t1 != null && t2 != null && !t1.started && !t2.started) {
                int index = taskQueue.indexOf(t2);
                taskQueue.remove(t1);
                taskQueue.add(index, t1);
            }
        }
    }

    public void moveTasksDown(ArrayList<Integer> taskIDs) {
        EllaTask t1 = null, t2 = null;
        for (int id : taskIDs) {
            for (EllaTask t : taskQueue) {
                if (t.getId() == id) {
                    t1 = t;
                } else if (t1 != null) {
                    t2 = t;
                    break;
                }
            }
            if (t1 != null && t2 != null && !t1.started && !t2.started) {
                int index = taskQueue.indexOf(t1);
                taskQueue.remove(t2);
                taskQueue.add(index, t2);
            }
        }
    }

    public class EllaTask implements Runnable {

        private SimpleStringProperty status = new SimpleStringProperty("waiting");
        private boolean started = false, canceled = false;
        private RuntimeThread runtimeThread;
        private ProfileThread profileThread;
        private int id;
        private String mode;
        private String settings, description;
        private DoubleProperty progress = new SimpleDoubleProperty(0.);
        private Aligner aligner;
        private IndexCreator indexCreator;

        public EllaTask(int id, String mode, String description, String settings) {
            this.description = description;
            this.id = id;
            this.runtimeThread = new RuntimeThread(id, presenter);
            this.profileThread = new ProfileThread(id, presenter);
            this.mode = mode;
            this.settings = settings;
        }

        public void run() {

            try {

                status.set("running");
                runtimeThread.start();
                profileThread.start();
                runningTask = this;

                String[] args = (mode + " " + settings).split("\\s+");
                switch (mode) {
                    case "blastx":
                        AlignOptionHandler.run(args, true);
                        File indexFile = AlignOptionHandler.indexFile;
                        File queryFile = AlignOptionHandler.queryFile;
                        int m = AlignOptionHandler.m;
                        int alignCores = AlignOptionHandler.cores;
                        File daaFile = AlignOptionHandler.outputFile;
                        aligner = new Aligner(indexFile, queryFile, daaFile, m, alignCores);
                        aligner.setTaskmanager(Taskmanager.this);
                        progress.bind(aligner.totalProgressProperty());
                        presenter.reportAligner(id, aligner);
                        aligner.run();
                        runningTaskProgress.bind(aligner.totalProgressProperty());
                        break;
                    case "index":
                        IndexOptionHandler.run(args);
                        File inFile = IndexOptionHandler.inFile;
                        File dbFile = IndexOptionHandler.dbFile;
                        int cores = IndexOptionHandler.cores;
                        Long size = IndexOptionHandler.size;
                        indexCreator = new IndexCreator();
                        indexCreator.setTaskmanager(Taskmanager.this);
                        progress.bind(indexCreator.totalProgressProperty());
                        presenter.reportIndexCreator(id, indexCreator);
                        indexCreator.run(inFile, dbFile, 3, MyParameters.NO_SEED_SHAPE, 3, 3, cores, size);
                        runningTaskProgress.bind(indexCreator.totalProgressProperty());
                        break;
                }

                if (!canceled) {
                    status.set("completed");
                    setStatus("completed");
                }

            } catch (Exception e) {
                e.printStackTrace();
            }
        }

        public void cancel() {
            status.set("canceled");
            setStatus("canceled");
            canceled = true;
            stop();
        }

        public void stop() {
            if (aligner != null)
                aligner.stop();
            if (indexCreator != null)
                indexCreator.stop();
            runtimeThread.setDone(true);
            profileThread.setDone(true);
        }

        public void start() {
            if (!started) {
                executor.submit(this);
                started = true;
            }
        }

        public int getId() {
            return id;
        }

        public String getMode() {
            return mode;
        }

        public String getSettings() {
            return settings;
        }

        public double getProgress() {
            return progress.get();
        }

        public DoubleProperty progressProperty() {
            return progress;
        }

        public String getStatus() {
            return status.get();
        }

        public SimpleStringProperty statusProperty() {
            return status;
        }

        public String getDescription() {
            return description;
        }
    }

    public ObservableList<EllaTask> getTaskQueue() {
        return taskQueue;
    }

    public boolean isIsRunning() {
        return isRunning.get();
    }

    public BooleanProperty isRunningProperty() {
        return isRunning;
    }

    public String getStatus() {
        return status.get();
    }

    public StringProperty statusProperty() {
        return status;
    }

    public void setPresenter(Presenter presenter) {
        this.presenter = presenter;
    }

    public DoubleProperty runningTaskProgressProperty() {
        return runningTaskProgress;
    }
}
