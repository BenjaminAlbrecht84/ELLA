package ella.view.taskView;

import ella.presenter.Presenter;
import ella.view.EllaView;
import ella.view.detailsView.DetailsView;
import javafx.beans.property.ReadOnlyIntegerProperty;
import javafx.collections.ObservableList;
import javafx.scene.Scene;
import javafx.scene.control.ContextMenu;
import javafx.scene.control.MenuItem;
import javafx.scene.input.MouseButton;
import javafx.scene.input.MouseEvent;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;

import java.util.HashMap;

public class TaskView extends BorderPane {

    private HashMap<Integer, Stage> id2detailsStage = new HashMap<>();
    private EllaView view;
    private TaskViewTable taskViewTable;

    private MenuItem removeItem, cancelItem;
    private ContextMenu contextMenu = new ContextMenu();
    private Presenter presenter;

    public TaskView(EllaView view) {
        this.view = view;
        initFields();
        setupLayout();
        setUpBindings(view);
    }

    private void setUpBindings(EllaView view) {

        prefHeightProperty().bind(view.getRoot().heightProperty());
        prefWidthProperty().bind(view.getRoot().widthProperty());

        taskViewTable.prefHeightProperty().bind(view.getRoot().heightProperty());
        taskViewTable.prefWidthProperty().bind(view.getRoot().widthProperty());

        taskViewTable.addEventHandler(MouseEvent.MOUSE_CLICKED, t -> {
            if (t.getButton() == MouseButton.SECONDARY)
                contextMenu.show(taskViewTable, t.getScreenX(), t.getScreenY());
        });
    }

    public void setPresenter(Presenter presenter) {
        this.presenter = presenter;
        taskViewTable.setTaskmanager(presenter.getTaskmanager());
    }

    private void initFields() {
        taskViewTable = new TaskViewTable(this);
    }

    private void setupLayout() {

        // adding task table
        setCenter(taskViewTable);

        // setting up context menu
        removeItem = new MenuItem("Remove");
        cancelItem = new MenuItem("Cancel");
        contextMenu.getItems().addAll(cancelItem, removeItem);

    }

    public void addLoads(int id, String time, double cpuLoad, double usedMemory, double allocatedMemory) {
        DetailsView detailsView = (DetailsView) id2detailsStage.get(id).getScene().getRoot();
        detailsView.addLoads(time, cpuLoad, usedMemory, allocatedMemory);
    }

    public void showDetailsView(int id) {
        id2detailsStage.get(id).show();
    }

    public void updateRuntime(int id, String runtime) {
        DetailsView detailsView = (DetailsView) id2detailsStage.get(id).getScene().getRoot();
        detailsView.updateRuntime(runtime);
    }

    public void reportLogInfo(Integer id, String info) {
        DetailsView detailsView = (DetailsView) id2detailsStage.get(id).getScene().getRoot();
        detailsView.reportLogInfo(info);
    }

    public void reportNewTask(int id, String mode, String settings) {

        // setting up details view
        Stage detailsStage = new Stage();
        detailsStage.setScene(new Scene(new DetailsView()));
        id2detailsStage.put(id, detailsStage);
        DetailsView detailsView = (DetailsView) id2detailsStage.get(id).getScene().getRoot();
        detailsView.reportSettings(mode, settings);

    }

    public void reportProgressList(int id, ObservableList<Object[]> progressList) {
        DetailsView detailsView = (DetailsView) id2detailsStage.get(id).getScene().getRoot();
        detailsView.reportProgressList(progressList);
    }

    public ReadOnlyIntegerProperty selectedTaskProperty() {
        return taskViewTable.getSelectionModel().selectedIndexProperty();
    }

    public void deleteSelectedTasks() {
        presenter.removeTask(taskViewTable.getSelectedTaskIDs());
    }

    public void cancelSelectedTasks() {
        presenter.cancelTasks(taskViewTable.getSelectedTaskIDs());
    }

    public void selectAllTasks(boolean b) {
        if (b)
            taskViewTable.getSelectionModel().selectAll();
        else
            taskViewTable.getSelectionModel().clearSelection();
    }

    public void moveTask(boolean moveUpwards) {
        presenter.moveTaks(taskViewTable.getSelectedTaskIDs(), moveUpwards);
    }

}
