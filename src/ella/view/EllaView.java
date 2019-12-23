package ella.view;

import ella.view.convertView.ConvertController;
import ella.view.convertView.ConvertView;
import ella.view.settings.index_settings.IndexSettings;
import ella.view.settings.index_settings.IndexSettingsController;
import ella.view.taskView.TaskView;
import ella.view.utils.IconFactory;
import javafx.beans.property.DoubleProperty;
import javafx.collections.ObservableList;
import javafx.fxml.FXMLLoader;
import javafx.geometry.Pos;
import javafx.geometry.Rectangle2D;
import javafx.scene.Scene;
import javafx.scene.control.*;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;
import javafx.scene.text.TextAlignment;
import javafx.stage.Screen;
import javafx.stage.Stage;
import ella.presenter.Presenter;

import java.text.DecimalFormat;
import java.time.Duration;

import ella.view.settings.align_settings.AlignSettings;
import ella.view.settings.align_settings.AlignSettingsController;

/**
 * Created by Benjamin on 28.09.18.
 */
public class EllaView {

    private BorderPane root;
    private Presenter presenter;

    private DecimalFormat decimalFormat = new DecimalFormat("#.00");

    private ProgressBar progress;
    private Label status, computation;

    private Pane toolbar;
    private Button upButton = new Button(), downButton = new Button(), runButton = new Button(), pauseButton = new Button(), cancelButton = new Button(), deleteButton = new Button(), addButton = new Button(), convertButton = new Button();

    private TaskView taskView;

    private Stage alignStage, indexStage, convertStage;
    private AlignSettings alignSettings;
    private IndexSettings indexSettings;
    private ConvertView convertView;

    private MenuItem alignItem, indexItem, closeItem;
    private MenuItem runItem, pauseItem;
    private MenuItem convertItem;
    private MenuItem cancelItem, deleteItem, selectItem, deselectItem, upItem, downItem;
    private MenuItem aboutItem;

    public EllaView() {
        initFields();
        setUpLayout();
        setUpActions();
        setUpBindings();
    }

    private void setUpBindings() {
        cancelItem.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        cancelButton.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        deleteItem.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        deleteButton.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        upItem.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        upButton.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        downItem.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
        downButton.disableProperty().bind(taskView.selectedTaskProperty().isEqualTo(-1));
    }

    private void setUpActions() {
        alignItem.setOnAction(e -> alignSettings.show());
        addButton.setOnAction(e -> alignSettings.show());
        indexItem.setOnAction(e -> indexSettings.show());
        convertItem.setOnAction(e -> convertView.show());
        convertButton.setOnAction(e -> convertView.show());
        runItem.setOnAction(e -> presenter.pauseResumeTasks());
        runButton.setOnAction(e -> presenter.pauseResumeTasks());
        pauseItem.setOnAction(e -> presenter.pauseResumeTasks());
        pauseButton.setOnAction(e -> presenter.pauseResumeTasks());
        cancelItem.setOnAction(e -> taskView.cancelSelectedTasks());
        cancelButton.setOnAction(e -> taskView.cancelSelectedTasks());
        selectItem.setOnAction(e -> taskView.selectAllTasks(true));
        deselectItem.setOnAction(e -> taskView.selectAllTasks(false));
        deleteItem.setOnAction(e -> taskView.deleteSelectedTasks());
        deleteButton.setOnAction(e -> taskView.deleteSelectedTasks());
        upItem.setOnAction(e -> taskView.moveTask(true));
        upButton.setOnAction(e -> taskView.moveTask(true));
        downItem.setOnAction(e -> taskView.moveTask(false));
        downButton.setOnAction(e -> taskView.moveTask(false));
        aboutItem.setOnAction(e -> showSplashScreen());

    }

    private void setUpLayout() {

        Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();
        root.setPrefSize(primaryScreenBounds.getWidth() * 0.75, primaryScreenBounds.getHeight() * 0.75);
        root.setMinSize(primaryScreenBounds.getWidth() * 0.75, primaryScreenBounds.getHeight() * 0.75);

        // setting task tab
        taskView = new TaskView(this);
        root.setCenter(taskView);

        status.prefWidthProperty().bind(root.widthProperty().divide(3));
        progress.prefWidthProperty().bind(root.widthProperty().divide(3).add(-50));

        computation.prefWidthProperty().bind(root.widthProperty().divide(3));
        computation.setTextAlignment(TextAlignment.RIGHT);
        computation.setAlignment(Pos.CENTER_RIGHT);
        computation.setContentDisplay(ContentDisplay.RIGHT);

    }

    private void initFields() {

        try {

            // initializing main window
            FXMLLoader ellaViewloader = new FXMLLoader(EllaController.class.getResource("EllaWindow.fxml"));
            ellaViewloader.load();
            EllaController controller = ellaViewloader.getController();

            root = ellaViewloader.getRoot();

            // initializing file menu items
            this.alignItem = controller.alignItem;
            this.indexItem = controller.indexItem;
            this.closeItem = controller.closeItem;
            this.runItem = controller.runItem;
            this.pauseItem = controller.pauseItem;

            // initializing view menu items
            this.convertItem = controller.convertItem;

            // initializing edit menu items
            this.cancelItem = controller.cancelItem;
            this.selectItem = controller.selectItem;
            this.deleteItem = controller.delItem;
            this.deselectItem = controller.deselectItem;
            this.upItem = controller.upItem;
            this.downItem = controller.downItem;

            // initializing help menu items
            this.aboutItem = controller.aboutItem;

            // initializing toolbar
            this.toolbar = controller.toolbar;
            HBox toolbox = new HBox();
            toolbox.setSpacing(1);
            addButton.setGraphic(IconFactory.createImageView("icons8-add-file-50.png", toolbar.heightProperty()));
            Tooltip.install(addButton, new Tooltip("Add new align task."));
            runButton.setGraphic(IconFactory.createImageView("icons8-start-50.png", toolbar.heightProperty()));
            Tooltip.install(runButton, new Tooltip("Run all ELLA tasks."));
            pauseButton.setGraphic(IconFactory.createImageView("icons8-pause-50.png", toolbar.heightProperty()));
            Tooltip.install(pauseButton, new Tooltip("Pause all ELLA tasks."));
            cancelButton.setGraphic(IconFactory.createImageView("icons8-cancel-50.png", toolbar.heightProperty()));
            Tooltip.install(cancelButton, new Tooltip("Cancel a running ELLA task."));
            deleteButton.setGraphic(IconFactory.createImageView("icons8-delete-trash-50.png", toolbar.heightProperty()));
            Tooltip.install(deleteButton, new Tooltip("Delete an ELLA task."));
            convertButton.setGraphic(IconFactory.createImageView("icons8-rich-text-converter-50.png", toolbar.heightProperty()));
            Tooltip.install(convertButton, new Tooltip("Convert a DAA file into BLAST format."));
            upButton.setGraphic(IconFactory.createImageView("icons8-scroll-up-48.png", toolbar.heightProperty()));
            Tooltip.install(upButton, new Tooltip("Move selected ELLA task upwards."));
            downButton.setGraphic(IconFactory.createImageView("icons8-below-50.png", toolbar.heightProperty()));
            Tooltip.install(downButton, new Tooltip("Move selected ELLA task downwards."));
            toolbox.getChildren().addAll(addButton, new Separator(), runButton, pauseButton, new Separator(), upButton, downButton, cancelButton, deleteButton, new Separator(), convertButton);
            toolbar.getChildren().add(toolbox);

            // initializing bottom pane
            this.status = controller.status;
            this.progress = controller.progress;
            this.computation = controller.computation;

            //initializing align settings window
            FXMLLoader alignSettingsLoader = new FXMLLoader(AlignSettingsController.class.getResource("AlignSettingsWindow.fxml"));
            alignSettingsLoader.load();
            alignStage = new Stage();
            alignStage.setTitle("ELLAs Alignment Setup");
            alignSettings = new AlignSettings(alignStage, alignSettingsLoader.getController());
            alignStage.setScene(new Scene(alignSettingsLoader.getRoot()));

            //initializing index settings window
            FXMLLoader indexSettingsLoader = new FXMLLoader(IndexSettingsController.class.getResource("IndexSettingsWindow.fxml"));
            indexSettingsLoader.load();
            indexStage = new Stage();
            indexStage.setTitle("ELLAs Index Setup");
            indexSettings = new IndexSettings(indexStage, indexSettingsLoader.getController());
            indexStage.setScene(new Scene(indexSettingsLoader.getRoot()));

            // initializing convert window
            FXMLLoader convertLoader = new FXMLLoader(ConvertController.class.getResource("ConvertWindow.fxml"));
            convertLoader.load();
            convertStage = new Stage();
            convertStage.setTitle("ELLAs DAA Converter");
            convertView = new ConvertView(convertStage, convertLoader.getController());
            convertStage.setScene(new Scene(convertLoader.getRoot()));

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void addLoads(int id, String time, double cpuLoad, double usedMemory, double allocatedMemory) {
        computation.setText("CPU: " + (int) cpuLoad + "%   Memory: " + decimalFormat.format(usedMemory) + "GB");
        taskView.addLoads(id, time, cpuLoad, usedMemory, allocatedMemory);
    }

    public void reportLogInfo(Integer id, String info) {
        taskView.reportLogInfo(id, info);
    }

    public void reportNewTask(int id, String mode, String settings) {
        taskView.reportNewTask(id, mode, settings);
    }

    public void reportProgressList(int id, ObservableList<Object[]> progressList) {
        taskView.reportProgressList(id, progressList);
    }

    public void finish() {
        computation.setText("CPU: -   Memory: -");
    }

    public void setPresenter(Presenter presenter) {
        this.presenter = presenter;
        taskView.setPresenter(presenter);
        alignSettings.setPresenter(presenter);
        indexSettings.setPresenter(presenter);
        convertView.setPresenter(presenter);
        runItem.disableProperty().bind(presenter.isRunningProperty());
        runButton.disableProperty().bind(presenter.isRunningProperty());
        pauseItem.disableProperty().bind(presenter.isRunningProperty().not());
        pauseButton.disableProperty().bind(presenter.isRunningProperty().not());
    }

    public void updateRuntime(int id, String runtime) {
        taskView.updateRuntime(id, runtime);
    }

    public BorderPane getRoot() {
        return root;
    }

    public MenuItem getCloseItem() {
        return closeItem;
    }

    public Label getStatus() {
        return status;
    }

    public ProgressBar getProgress() {
        return progress;
    }

    public void bindProgress(DoubleProperty totalProgressProperty) {
        progress.progressProperty().bind(totalProgressProperty.divide(100.));
    }

    public void reportConvertStatus(String message, boolean isRunning) {
        convertView.reportConvertStatus(message, isRunning);
    }

    public void showSplashScreen() {

    }

}
