package ella.presenter;

import ella.model.Taskmanager;
import ella.model.aligner.aligning.Aligner;
import ella.model.aligner.indexing.creator1.IndexCreator;
import ella.model.io.converting.DaaConverter;
import javafx.application.Platform;
import ella.view.EllaView;
import javafx.beans.property.BooleanProperty;

import java.util.ArrayList;

/**
 * Created by Benjamin on 28.09.18.
 */
public class Presenter {

    private Taskmanager taskmanager;
    private EllaView view;
    private Thread convertThread;

    public Presenter(Taskmanager taskmanager, EllaView view) {
        this.taskmanager = taskmanager;
        this.view = view;
        setOnAction();
        initBindings(view);
    }

    private void initBindings(EllaView view) {
        view.getStatus().textProperty().bind(taskmanager.statusProperty());
        view.getProgress().progressProperty().bind(taskmanager.runningTaskProgressProperty());
    }

    public void addEllaTask(String settings, String description, String mode) {
        taskmanager.addTask(settings, description, mode);
    }

    public void pauseResumeTasks() {
        taskmanager.pauseResume();
    }

    public void cancelTasks(ArrayList<Integer> ids) {
        taskmanager.cancelTasks(ids);
    }

    public void removeTask(ArrayList<Integer> ids) {
        taskmanager.removeTasks(ids);
    }

    private void setOnAction() {
        view.getCloseItem().setOnAction(e -> {
            System.exit(0);
        });
    }

    public void reportAligner(int id, Aligner aligner) {
        view.reportProgressList(id, aligner.getProgressList());
        view.bindProgress(aligner.totalProgressProperty());
    }

    public void reportIndexCreator(int id, IndexCreator indexCreator) {
        view.bindProgress(indexCreator.totalProgressProperty());
    }

    public void reportRuntime(int id, String runtime) {
        Platform.runLater(() -> view.updateRuntime(id, runtime));
    }

    public void reportCpuLoad(int id, String time, double cpuLoad, double usedMemory, double allocatedMemory) {
        Platform.runLater(() -> view.addLoads(id, time, cpuLoad, usedMemory, allocatedMemory));
    }

    public void finish() {
        view.finish();
    }

    public Taskmanager getTaskmanager() {
        return taskmanager;
    }

    public void reportLogInfo(Integer id, String info) {
        view.reportLogInfo(id, info);
    }

    public void reportNewTask(int id, String mode, String settings, String description) {
        view.reportNewTask(id, mode, settings, description);
    }

    public int getNextTaskID() {
        return taskmanager.getNextTaskId();
    }

    public void convertFile(String daaFile, String outFile, String format, String cores) {
        convertThread = new Thread(() -> new DaaConverter().run(this, daaFile, outFile, format, cores));
        convertThread.start();
    }

    public void cancelConvertTask() {
        if (convertThread != null)
            convertThread.interrupt();
    }

    public BooleanProperty isRunningProperty() {
        return taskmanager.isRunningProperty();
    }

    public void reportConvertStatus(String message, boolean isRunning) {
        view.reportConvertStatus(message, isRunning);
    }

    public void moveTaks(ArrayList<Integer> taskIDs, boolean moveUpwards) {
        if (moveUpwards)
            taskmanager.moveTasksUp(taskIDs);
        else
            taskmanager.moveTasksDown(taskIDs);
    }
}
