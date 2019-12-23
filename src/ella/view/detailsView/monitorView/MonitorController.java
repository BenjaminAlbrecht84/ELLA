package ella.view.detailsView.monitorView;

import javafx.fxml.FXML;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.NumberAxis;
import javafx.scene.chart.StackedAreaChart;
import javafx.scene.control.Label;
import javafx.scene.layout.HBox;

public class MonitorController {

    @FXML
    protected Label uptime;

    @FXML
    protected Label workLabel;

    @FXML
    protected HBox workflow;

    @FXML
    protected LineChart cpu;

    @FXML
    protected HBox cpuBox;

    @FXML
    protected HBox heapBox;

    @FXML
    protected NumberAxis xCPU;

    @FXML
    protected StackedAreaChart heap;

}
