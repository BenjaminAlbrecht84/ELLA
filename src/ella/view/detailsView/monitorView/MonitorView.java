package ella.view.detailsView.monitorView;

import ella.view.EllaController;
import ella.view.EllaView;
import ella.view.detailsView.DetailsView;
import javafx.application.Platform;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.ObservableList;
import javafx.fxml.FXMLLoader;
import javafx.scene.chart.LineChart;
import javafx.scene.chart.StackedAreaChart;
import javafx.scene.chart.XYChart;
import javafx.scene.control.Label;
import javafx.scene.layout.BorderPane;
import javafx.scene.layout.HBox;

import javax.management.monitor.Monitor;
import java.text.DecimalFormat;

public class MonitorView {

    private DecimalFormat decimalFormat = new DecimalFormat("#.00");

    private BorderPane root;
    private MyWorkflow workflow;
    private Label workLabel, uptime;

    private LineChart cpu;
    private Label cpuNum = new Label(), cpuUsage = new Label();
    private SimpleStringProperty cpuNumInfo = new SimpleStringProperty(), cpuUsageInfo = new SimpleStringProperty();
    private XYChart.Series<String, Double> cpuSeries = new XYChart.Series<String, Double>();
    private HBox cpuBox;

    private SimpleStringProperty usedMemInfo = new SimpleStringProperty(), allocMemInfo = new SimpleStringProperty();
    private Label usedMem = new Label(), allocMem = new Label();
    private StackedAreaChart heap;
    private HBox heapBox;
    private XYChart.Series<String, Double> usedMemSeries = new XYChart.Series<String, Double>();
    private XYChart.Series<String, Double> allocMemSeries = new XYChart.Series<String, Double>();

    public MonitorView(DetailsView detailsView) {
        initFields();
        root.prefHeightProperty().bind(detailsView.heightProperty().add(-25));
        root.prefWidthProperty().bind(detailsView.widthProperty());
    }

    private void initFields() {

        try {

            FXMLLoader loader = new FXMLLoader(MonitorController.class.getResource("MonitorView.fxml"));
            loader.load();
            root = loader.getRoot();

            MonitorController controller = loader.getController();

            this.uptime = controller.uptime;
            this.workflow = new MyWorkflow(controller.workflow);
            this.workLabel = controller.workLabel;
            this.cpu = controller.cpu;
            this.cpuBox = controller.cpuBox;
            cpu.getData().add(cpuSeries);
            this.heap = controller.heap;
            this.heapBox = controller.heapBox;
            heap.getData().addAll(usedMemSeries, allocMemSeries);

            cpuBox.getChildren().add(cpuUsage);
            cpuUsage.textProperty().bind(cpuUsageInfo);

            heapBox.getChildren().addAll(usedMem, allocMem);
            usedMem.textProperty().bind(usedMemInfo);
            allocMem.textProperty().bind(allocMemInfo);

        } catch (Exception e) {
            e.printStackTrace();
        }

    }

    public void addLoads(String time, double cpuLoad, double usedMemory, double allocatedMemory) {
        cpuSeries.getData().add(new XYChart.Data<>(time, cpuLoad));
        usedMemSeries.getData().add(new XYChart.Data<>(time, usedMemory));
        double cpuAverage = cpuSeries.getData().stream().mapToDouble(XYChart.Data::getYValue).average().getAsDouble();
        double usedMemMax = usedMemSeries.getData().stream().mapToDouble(XYChart.Data::getYValue).max().getAsDouble();
        cpuUsageInfo.set("Usage: " + (int) cpuLoad + "% (avg " + decimalFormat.format(cpuAverage) + "%)");
        usedMemInfo.set("Used: " + decimalFormat.format(usedMemory) + "GB (max " + decimalFormat.format(usedMemMax) + "GB)");
        allocMemInfo.set("Size: " + decimalFormat.format(allocatedMemory) + "GB ");
        allocMemSeries.getData().add(new XYChart.Data<>(time, allocatedMemory - usedMemory));
    }

    public BorderPane getRoot() {
        return root;
    }

    public void setCpuNum(int cpus) {
        cpuNum.setText(String.valueOf(cpus));
    }

    public void updateRuntime(String runningTime) {
        uptime.setText("Uptime : " + runningTime);
    }

    public void setUpWorkflow(ObservableList<Object[]> progressList) {
        workflow.setUp(progressList);
    }


}
