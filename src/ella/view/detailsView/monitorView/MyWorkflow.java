package ella.view.detailsView.monitorView;

import javafx.application.Platform;
import javafx.beans.property.SimpleDoubleProperty;
import javafx.collections.ListChangeListener;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.geometry.Pos;
import javafx.scene.control.ProgressBar;
import javafx.scene.control.ProgressIndicator;
import javafx.scene.control.TitledPane;
import javafx.scene.layout.*;

public class MyWorkflow {

    private HBox queue;

    public MyWorkflow(HBox queue) {
        this.queue = queue;
        queue.setPadding(new Insets(10, 10, 10, 10));
    }

    public void setUp(ObservableList<Object[]> data) {
        data.addListener((ListChangeListener<Object[]>) c -> {
            while (c.next()) {
                if (c.wasAdded()) {
                    for (Object[] o : c.getAddedSubList()) {
                        String title = (String) o[0];
                        String name1 = (String) o[1];
                        SimpleDoubleProperty progress1 = (SimpleDoubleProperty) o[2];
                        String name2 = (String) o[3];
                        SimpleDoubleProperty progress2 = (SimpleDoubleProperty) o[4];
                        addProgreesPane(title, name1, progress1, name2, progress2);
                    }
                }
            }
        });
    }

    private void addProgreesPane(String title, String name1, SimpleDoubleProperty progress1, String name2, SimpleDoubleProperty progress2) {
        VBox vBox = new VBox(addProgressNode(name1, progress1), addProgressNode(name2, progress2));
        vBox.setPadding(new Insets(10, 10, 10, 10));
        TitledPane pane = new TitledPane(title, vBox);
        pane.setCollapsible(false);
        Platform.runLater(() -> queue.getChildren().add(pane));
    }

    private TitledPane addProgressNode(String name, SimpleDoubleProperty progress) {
        ProgressBar progressBar = new ProgressBar();
        progressBar.progressProperty().bind(progress);
        ProgressIndicator progressIndicator = new ProgressIndicator();
        progressIndicator.progressProperty().bind(progress);
        HBox hBox = new HBox(progressBar, progressIndicator);
        hBox.setPadding(new Insets(10, 10, 10, 10));
        hBox.setSpacing(5);
        hBox.setAlignment(Pos.CENTER);
        TitledPane pane = new TitledPane(name, hBox);
        pane.setCollapsible(false);
        return pane;
    }

}
