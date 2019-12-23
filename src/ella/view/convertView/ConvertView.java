package ella.view.convertView;

import ella.model.io.converting.DaaConverter;
import ella.presenter.Presenter;
import javafx.application.Platform;
import javafx.beans.property.StringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.Node;
import javafx.scene.control.*;
import javafx.scene.layout.Priority;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

import java.io.*;
import java.util.List;
import java.util.Properties;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

public class ConvertView {

    private Pattern coresPattern = Pattern.compile("\\d*|\\d+");
    private TextFormatter coresFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> coresPattern.matcher(change.getControlNewText()).matches() ? (change.getText().isEmpty() || Integer.parseInt(change.getText()) > 0 ? change : null) : null);

    private Properties properties;
    private TextField coresField;
    private TextArea daaArea, outArea;
    private Button convertButton, outButton, daaButton;
    private Label infoLabel;

    private ComboBox<String> formatBox;
    private ObservableList<String> formatList = FXCollections.observableArrayList(DaaConverter.BLAST_FORMATS);

    private Presenter presenter;
    private Stage stage;

    public ConvertView(Stage stage, ConvertController controller) {

        this.stage = stage;
        this.properties = readInProperties();

        daaArea = controller.daaArea;
        outArea = controller.outArea;
        coresField = controller.coresField;

        daaButton = controller.daaButton;
        outButton = controller.outButton;
        convertButton = controller.convertButton;

        infoLabel = controller.infoLabel;

        formatBox = controller.formatBox;

        setUpActions();
        setUpLayout();
    }

    private void setUpActions() {
        daaButton.setOnAction(e -> {
            chooseFiles(daaArea.textProperty(), "Select DAA file", "DaaFile");
            setOutputFiles();
        });
        outButton.setOnAction(e -> chooseFolder(outArea.textProperty(), "Select output folder"));
        convertButton.setOnAction(e -> {
            String[] daas = daaArea.getText().split("\n");
            String[] outputs = outArea.getText().split("\n");
            if (daas.length == outputs.length) {
                for (int i = 0; i < daas.length; i++) {
                    String d = new File(daas[i]).getName();
                    String o = new File(outputs[i]).getName();

                    presenter.convertFile(daas[i], outputs[i], formatBox.getSelectionModel().getSelectedItem(), coresField.getText());
                    convertButton.setDisable(true);
                    infoLabel.setText("Converting " + d + " to " + o + "...");
                }
            }
        });
    }

    private void setUpLayout() {
        formatBox.setItems(formatList);
        formatBox.getSelectionModel().selectFirst();
        daaArea.setText(properties.getProperty("DaaFile"));
        if (!daaArea.getText().isEmpty())
            setOutputFiles();
        coresField.setText("1");
        coresField.setTextFormatter(coresFormatter);
    }

    private void setOutputFiles() {
        String[] daas = daaArea.getText().split("\n");
        outArea.setText("");
        for (String daa : daas) {
            String ending = formatBox.getSelectionModel().getSelectedItem().split("\\*.")[1].replaceAll("\\)", "");
            outArea.appendText(daa.replace(getFileEnding(daa), "." + ending));
        }
    }

    private String getFileEnding(String path) {
        int lastIndexOf = path.lastIndexOf(".");
        if (lastIndexOf == -1)
            return "";
        return path.substring(lastIndexOf);
    }

    private void chooseFiles(StringProperty prop, String title, String key) {
        FileChooser chooser = new FileChooser();
        chooser.setTitle(title);
        List<File> dir = chooser.showOpenMultipleDialog(stage);
        if (dir != null) {
            StringBuilder s = new StringBuilder();
            for (File f : dir)
                s.append(f.getAbsolutePath() + "\n");
            prop.set(s.toString());
            setProperties(key, s.toString());
        }
    }

    private void chooseFolder(StringProperty s, String title) {
        DirectoryChooser chooser = new DirectoryChooser();
        chooser.setTitle(title);
        File dir = chooser.showDialog(stage);
        if (dir != null) {
            s.set(dir.getAbsolutePath());
            setProperties(title, dir.getAbsolutePath());
        }
    }

    private void setProperties(String key, String value) {
        properties.setProperty(key, value);
        writeOutProperties();
    }

    private Properties readInProperties() {
        properties = new Properties();
        try {
            BufferedInputStream stream = new BufferedInputStream(new FileInputStream("ella.properties"));
            properties.load(stream);
            stream.close();
        } catch (Exception e) {
            writeOutProperties();
        }
        return properties;
    }

    private void writeOutProperties() {
        try {
            BufferedOutputStream stream = new BufferedOutputStream(new FileOutputStream("ella.properties"));
            properties.store(stream, "ELLAs specific properties");
            stream.close();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

    public void setPresenter(Presenter presenter) {
        this.presenter = presenter;
        stage.setOnCloseRequest(e -> presenter.cancelConvertTask());
    }

    public void show() {
        infoLabel.setText("");
        convertButton.setDisable(false);
        stage.show();
    }

    public void reportConvertStatus(String message, boolean isRunning) {
        if (!isRunning)
            convertButton.setDisable(false);
        Platform.runLater(() -> infoLabel.setText(message));
    }

}
