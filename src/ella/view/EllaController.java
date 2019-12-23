package ella.view;

import javafx.fxml.FXML;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedAreaChart;
import javafx.scene.control.*;
import javafx.scene.layout.HBox;
import javafx.scene.layout.Pane;


/**
 * Created by Benjamin on 01.10.18.
 */
public class EllaController {

    // File items **********************

    @FXML
    protected MenuItem alignItem, indexItem, runItem, pauseItem, closeItem;

    // Edit items **********************

    @FXML
    protected MenuItem selectItem, deselectItem, cancelItem, delItem, upItem, downItem;

    // View items **********************

    @FXML
    protected MenuItem convertItem;

    // Help items **********************

    protected MenuItem aboutItem;

    // **********************

    @FXML
    protected Pane toolbar;

    @FXML
    protected ProgressBar progress;

    @FXML
    protected Label computation;

    @FXML
    protected Label status;

}
