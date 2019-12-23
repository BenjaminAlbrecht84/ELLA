package ella.view.settings.index_settings;

import ella.presenter.Presenter;
import javafx.beans.property.StringProperty;
import javafx.scene.control.Button;
import javafx.scene.control.TextArea;
import javafx.scene.control.TextField;
import javafx.scene.control.TextFormatter;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;

import java.io.*;
import java.util.List;
import java.util.Properties;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

public class IndexSettings {

    private Pattern intPattern = Pattern.compile("\\d*|\\d+");
    private TextFormatter volumeFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> intPattern.matcher(change.getControlNewText()).matches() ? change : null);
    private TextFormatter coresFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> intPattern.matcher(change.getControlNewText()).matches() ? change : null);


    private Stage stage;
    private Properties properties;
    private Presenter presenter;

    private Button addButton, cancelButton;

    private TextArea input, database;
    private TextField cores, volume;
    private Button inButton, outButton;

    public IndexSettings(Stage stage, IndexSettingsController controller) {

        this.stage = stage;
        this.properties = readInProperties();

        this.addButton = controller.runButton;
        this.cancelButton = controller.cancelButton;

        this.input = controller.input;
        input.setEditable(false);
        this.inButton = controller.inButton;
        this.database = controller.output;
        this.outButton = controller.outButton;
        this.cores = controller.cores;
        this.volume = controller.volume;

        setUpActions();
        setUpBindings();
        setUpLayout();

    }

    private void setUpBindings() {
        outButton.disableProperty().bind(input.textProperty().isEmpty());
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

    private void setUpActions() {
        cancelButton.setOnAction(e -> stage.hide());
        inButton.setOnAction(e -> {
            chooseFiles(input.textProperty(), "Select reference protein file", "RefFile");
            setOutputFiles();
        });
        outButton.setOnAction(e -> {
            chooseOutputFolder(database.textProperty(), "Select database folder");
        });
        addButton.setOnAction(e -> {
            String[] inputs = input.getText().split("\n");
            String[] databases = database.getText().split("\n");
            if (inputs.length == databases.length) {
                for (int i = 0; i < inputs.length; i++) {
                    String in = new File(inputs[i]).getName();
                    String db = new File(databases[i]).getName();
                    presenter.addEllaTask(parseSettings(inputs[i], databases[i]), "indexing " + in + " into " + db, "index");
                }
                stage.hide();
            }
        });

    }

    private void setUpLayout() {
        input.setText(properties.getProperty("RefFile"));
        if (input.getText() != null && !input.getText().isEmpty())
            setOutputFiles();
        cores.setText(String.valueOf(Runtime.getRuntime().availableProcessors()));
        cores.setTextFormatter(coresFormatter);
        volume.setTextFormatter(volumeFormatter);
        volume.setText("");
    }

    private String parseSettings(String in, String db) {
        StringBuilder settings = new StringBuilder();
        settings.append("-i " + in);
        settings.append(" -d " + db);
        settings.append(" -p " + cores.getText());
        settings.append(" -s " + volume.getText() + "g");
        return settings.toString();
    }

    private void setOutputFiles() {
        String[] inputs = input.getText().split("\n");
        database.setText("");
        for (String in : inputs) {
            database.appendText(in.replace(getFileExtension(in), ".edb") + "\n");
        }
    }

    private void setOutputName(String input) {
        if (!input.isEmpty()) {
            String name = new File(input).getName();
            name = name.replace(getFileExtension(name), ".edb");
            database.setText(database.getText().concat(File.separatorChar + name));
        }
    }

    private String getFileExtension(String path) {
        int lastIndexOf = path.lastIndexOf(".");
        if (lastIndexOf == -1)
            return "";
        return path.substring(lastIndexOf);
    }

    private void chooseFile(StringProperty s, String title, String key) {
        FileChooser chooser = new FileChooser();
        chooser.setTitle(title);
        File dir = chooser.showOpenDialog(stage);
        if (dir != null) {
            s.set(dir.getAbsolutePath());
            setProperties(key, dir.getAbsolutePath());
        }
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

    private void chooseOutputFolder(StringProperty s, String title) {
        DirectoryChooser chooser = new DirectoryChooser();
        chooser.setTitle(title);
        File dir = chooser.showDialog(stage);
        if (dir != null && !input.getText().isEmpty()) {
            String[] files = input.getText().split("\n");
            StringBuilder content = new StringBuilder();
            for (String f : files) {
                if (!f.isEmpty()) {
                    String ending = getFileExtension(f);
                    String fileName = new File(f).getName().replaceAll(ending, ".edb");
                    content.append(dir.getAbsolutePath() + File.separator + fileName + "\n");
                }
            }
            s.setValue(content.toString());
        }
    }

    private void setProperties(String key, String value) {
        properties.setProperty(key, value);
        writeOutProperties();
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
    }

    public void show() {
        stage.show();
    }


}
