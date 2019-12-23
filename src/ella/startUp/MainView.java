package ella.startUp;

import ella.model.Taskmanager;
import javafx.application.Application;
import javafx.scene.Scene;
import javafx.scene.layout.BorderPane;
import javafx.stage.Stage;
import ella.presenter.Presenter;
import ella.view.EllaView;

import java.util.Locale;

/**
 * Created by Benjamin on 28.09.18.
 */
public class MainView extends Application {

    public static void main(String[] args) {
        Application.launch(args);
    }

    @Override
    public void start(Stage primaryStage) throws Exception {

        Locale.setDefault(new Locale("en", "US"));

        Taskmanager taskmanager = new Taskmanager();
        EllaView view = new EllaView();
        Presenter presenter = new Presenter(taskmanager, view);
        taskmanager.setPresenter(presenter);
        view.setPresenter(presenter);

        Scene scene = new Scene((BorderPane) view.getRoot());

        primaryStage.setScene(scene);
        primaryStage.setTitle("ELLA - Enhanced Local Alignment Tool");
        primaryStage.show();

        primaryStage.setOnCloseRequest(e -> {
            System.exit(0);
        });

    }

}
