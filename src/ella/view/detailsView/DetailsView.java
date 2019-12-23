package ella.view.detailsView;

import ella.view.detailsView.monitorView.MonitorView;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.geometry.Rectangle2D;
import javafx.scene.control.Tab;
import javafx.scene.control.TabPane;
import javafx.scene.control.TextArea;
import javafx.stage.Screen;

public class DetailsView extends TabPane {

    private SettingsView settingsView;
    private TextArea messageArea;
    private MonitorView monitorView;
    private Tab settingsTab, messageTab, monitorTab;

    public DetailsView() {
        initLayout();
    }

    private void initLayout() {

        Rectangle2D primaryScreenBounds = Screen.getPrimary().getVisualBounds();
        setPrefSize(primaryScreenBounds.getWidth() * 0.7, primaryScreenBounds.getHeight() * 0.75);

        // setting up monitor tab
        monitorTab = new Tab("Performance");
        monitorView = new MonitorView(this);
        monitorTab.setContent(monitorView.getRoot());

        // setting up settings tab
        settingsView = new SettingsView();
        settingsTab = new Tab("Settings");
        settingsTab.setContent(settingsView);

        // setting up log tab
        messageArea = new TextArea();
        messageTab = new Tab("Messages");
        messageTab.setContent(messageArea);

        getTabs().addAll(settingsTab, monitorTab, messageTab);
        setTabClosingPolicy(TabClosingPolicy.UNAVAILABLE);

    }

    public void addLoads(String time, double cpuLoad, double usedMemory, double allocatedMemory) {
        monitorView.addLoads(time, cpuLoad, usedMemory, allocatedMemory);
    }

    public void updateRuntime(String runtime) {
        monitorView.updateRuntime(runtime);
    }

    public void reportLogInfo(String info) {
        messageArea.appendText("\n" + info);
    }

    public void reportSettings(String mode, String settings) {
        settingsView.init(mode, settings, this);
    }

    public void reportProgressList(ObservableList<Object[]> progressList) {
        monitorView.setUpWorkflow(progressList);
    }

}
