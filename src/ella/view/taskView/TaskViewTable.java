package ella.view.taskView;

import javafx.scene.control.*;
import ella.model.Taskmanager;
import javafx.scene.control.cell.PropertyValueFactory;
import javafx.util.Callback;

import java.util.ArrayList;

public class TaskViewTable extends TableView<Taskmanager.EllaTask> {

    public TaskViewTable(TaskView taskView) {

        setPlaceholder(new Label("Please add some ELLA tasks: File -> New Align Task | New Index Task"));
        getSelectionModel().setSelectionMode(SelectionMode.MULTIPLE);

        prefWidthProperty().bind(taskView.widthProperty());
        prefHeightProperty().bind(taskView.heightProperty());
        setColumnResizePolicy(UNCONSTRAINED_RESIZE_POLICY);

        TableColumn<Taskmanager.EllaTask, Integer> idColumn = new TableColumn<>("Task ID");
        TableColumn<Taskmanager.EllaTask, String> modeColumn = new TableColumn<>("Mode");
        TableColumn<Taskmanager.EllaTask, String> descriptionColumn = new TableColumn<>("Description");
        TableColumn<Taskmanager.EllaTask, Double> statusColumn = new TableColumn<>("Status");
        TableColumn<Taskmanager.EllaTask, Double> progressColumn = new TableColumn<>("Progress");
        TableColumn detailsColumn = new TableColumn("");
        detailsColumn.setCellFactory(setUpButtonCellFactory(taskView, "Details"));
        getColumns().addAll(idColumn, modeColumn, descriptionColumn, statusColumn, progressColumn, detailsColumn);

        double n = getColumns().size() + 2;
        idColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(1. / n));
        idColumn.setStyle("-fx-alignment: CENTER;");
        modeColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(1. / n));
        modeColumn.setStyle("-fx-alignment: CENTER;");
        progressColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(1. / n));
        progressColumn.setStyle("-fx-alignment: CENTER;");
        detailsColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(1. / n));
        detailsColumn.setStyle("-fx-alignment: CENTER;");
        statusColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(1. / n));
        statusColumn.setStyle("-fx-alignment: CENTER;");
        descriptionColumn.prefWidthProperty().bind(taskView.widthProperty().multiply(3. / n));
        descriptionColumn.setStyle("-fx-alignment: CENTER;");

        statusColumn.setCellValueFactory(new PropertyValueFactory<>("status"));
        idColumn.setCellValueFactory(new PropertyValueFactory<>("id"));
        modeColumn.setCellValueFactory(new PropertyValueFactory<>("mode"));
        descriptionColumn.setCellValueFactory(new PropertyValueFactory<>("description"));
        progressColumn.setCellValueFactory(new PropertyValueFactory<>("progress"));

    }

    private Callback<TableColumn<Taskmanager.EllaTask, String>, TableCell<Taskmanager.EllaTask, String>> setUpButtonCellFactory(TaskView taskView, String text) {
        return new Callback<TableColumn<Taskmanager.EllaTask, String>, TableCell<Taskmanager.EllaTask, String>>() {
            @Override
            public TableCell call(final TableColumn<Taskmanager.EllaTask, String> param) {
                final TableCell<Taskmanager.EllaTask, String> cell = new TableCell<Taskmanager.EllaTask, String>() {

                    final Button button = new Button(text);

                    @Override
                    public void updateItem(String item, boolean empty) {
                        super.updateItem(item, empty);
                        if (empty) {
                            setGraphic(null);
                            setText(null);
                        } else {
                            Taskmanager.EllaTask task = getTableView().getItems().get(getIndex());
                            button.setOnAction(event -> taskView.showDetailsView(task.getId()));
                            setGraphic(button);
                            setText(null);
                        }
                    }
                };

                return cell;
            }
        };
    }

    public void setTaskmanager(Taskmanager taskmanager) {
        setItems(taskmanager.getTaskQueue());
    }

    public ArrayList<Integer> getSelectedTaskIDs() {
        ArrayList<Integer> taskIDs = new ArrayList<>();
        for (Taskmanager.EllaTask task : getSelectionModel().getSelectedItems()) {
            taskIDs.add(task.getId());
        }
        return taskIDs;
    }

}
