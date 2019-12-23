package ella.view.settings.align_settings;

import ella.presenter.Presenter;
import javafx.beans.property.StringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.scene.control.*;
import javafx.stage.DirectoryChooser;
import javafx.stage.FileChooser;
import javafx.stage.Stage;
import ella.model.aligner.utils.BlastStatisticsHelper;
import ella.model.io.MyParameters;

import java.io.*;
import java.util.List;
import java.util.Properties;
import java.util.function.UnaryOperator;
import java.util.regex.Pattern;

public class AlignSettings {

    private static String[] READ_TYPES = {"long reads", "short reads"};

    private Pattern doublePattern = Pattern.compile("\\d*|\\d+\\.\\d*");
    private TextFormatter evalueFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> doublePattern.matcher(change.getControlNewText()).matches() ? change : null);

    private Pattern intPattern = Pattern.compile("\\d*|\\d+");
    private TextFormatter multiplicityFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> intPattern.matcher(change.getControlNewText()).matches() ? change : null);
    private TextFormatter frameshiftFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> intPattern.matcher(change.getControlNewText()).matches() ? change : null);

    private Pattern coresPattern = Pattern.compile("\\d*|\\d+");
    private TextFormatter coresFormatter = new TextFormatter((UnaryOperator<TextFormatter.Change>) change -> coresPattern.matcher(change.getControlNewText()).matches() ? (change.getText().isEmpty() || Integer.parseInt(change.getText()) > 0 ? change : null) : null);

    private Stage stage;
    private Properties properties;
    private Presenter presenter;

    private Button addButton, cancelButton;
    private Label info;

    private TextArea input, output;
    private TextField reference, cores, multiplicity;
    private Button refButton, inButton, outButton;
    private ComboBox typeBox;

    private ComboBox matrices, penalties;
    private TextField frameshift;
    private ObservableList<String> matrixList = FXCollections.observableArrayList(BlastStatisticsHelper.matrixList);

    private Slider minCoverage;
    private Label minCovLabel;
    private TextField eValue;

    public AlignSettings(Stage stage, AlignSettingsController controller) {

        this.stage = stage;
        this.properties = readInProperties();

        this.addButton = controller.runButton;
        this.cancelButton = controller.cancelButton;
        this.info = controller.info;

        this.reference = controller.reference;
        reference.setEditable(false);
        this.refButton = controller.refButton;
        this.input = controller.inputArea;
        input.setEditable(false);
        this.inButton = controller.inButton;
        this.output = controller.outputArea;
        this.outButton = controller.outButton;
        this.cores = controller.cores;
        this.multiplicity = controller.multiplicity;
        this.typeBox = controller.typeBox;

        this.matrices = controller.matrices;
        this.penalties = controller.penalties;
        this.frameshift = controller.frameshift;

        this.minCoverage = controller.minCoverage;
        this.minCovLabel = controller.minCovLabel;
        this.eValue = controller.eValue;

        this.refButton = controller.refButton;

        setUpActions();
        setUpBindings();
        setUpLayout();

    }

    private void setUpBindings() {
        outButton.disableProperty().bind(input.textProperty().isEmpty());
    }

    private void setUpActions() {

        cancelButton.setOnAction(e -> stage.hide());
        matrices.setOnAction(e -> {
            adaptPenalties();
            updateMatrix();
        });
        refButton.setOnAction(e -> chooseFile(reference.textProperty(), "Select reference index file", "ReferenceFile"));
        inButton.setOnAction(e -> {
            chooseFiles(input.textProperty(), "Select input file", "InputFile");
            setOutputFiles();
        });
        outButton.setOnAction(e -> {
            chooseOutputFolder(output.textProperty(), "Select output folder");
        });
        addButton.setOnAction(e -> {
            String[] inputs = input.getText().split("\n");
            String[] outputs = output.getText().split("\n");
            if (inputs.length == outputs.length) {
                for (int i = 0; i < inputs.length; i++) {
                    String q = new File(inputs[i]).getName();
                    String r = new File(reference.getText()).getName();
                    String settings = parseSettings(inputs[i], outputs[i]);
                    presenter.addEllaTask(settings, "aligning " + q + " vs " + r, "blastx");
                }
                stage.hide();
            }
        });
        typeBox.valueProperty().addListener((observable, oldValue, newValue) -> {
            if (newValue.equals(READ_TYPES[0]))
                frameshift.setText("15");
            if (newValue.equals(READ_TYPES[1]))
                frameshift.setText("0");
        });

    }

    private String parseSettings(String in, String out) {
        StringBuilder settings = new StringBuilder();
        settings.append("-d " + reference.getText());
        settings.append(" -i " + in);
        settings.append(" -o " + out);
        settings.append(" -m " + multiplicity.getText());
        settings.append(" -p " + cores.getText());
        settings.append(" -F " + frameshift.getText());
        settings.append(" -e " + eValue.getText());
        settings.append(" -c " + (int) Math.round(minCoverage.getValue()));
        settings.append(" --matrix " + matrices.getSelectionModel().getSelectedItem());
        settings.append(" --gop " + ((String) penalties.getSelectionModel().getSelectedItem()).split("/")[0]);
        settings.append(" --gep " + ((String) penalties.getSelectionModel().getSelectedItem()).split("/")[1]);
        return settings.toString();
    }

    private void updateMatrix() {
        String type = (String) matrices.getSelectionModel().getSelectedItem();
        String[] split = ((String) penalties.getSelectionModel().getSelectedItem()).split("/");
        int gop = Integer.parseInt(split[0]);
        int gep = Integer.parseInt(split[1]);
        MyParameters.setScoringMatrix(type, gop, gep);
    }

    private void setOutputName(String input) {
        if (!input.isEmpty()) {
            String name = new File(input).getName();
            name = name.replace(getFileExtension(name), ".daa");
            output.setText(output.getText().concat(File.separatorChar + name));
        }
    }

    private void setOutputFiles() {
        String[] inputs = input.getText().split("\n");
        output.setText("");
        for (String in : inputs) {
            output.appendText(in.replace(getFileExtension(in), ".daa") + "\n");
        }
    }

    private String getFileExtension(String path) {
        int lastIndexOf = path.lastIndexOf(".");
        if (lastIndexOf == -1)
            return "";
        return path.substring(lastIndexOf);
    }


    private void adaptPenalties() {
        penalties.setItems(FXCollections.observableArrayList(BlastStatisticsHelper.matrix2penalties.get(matrices.getSelectionModel().getSelectedItem())));
        penalties.getSelectionModel().selectFirst();
    }

    private void setUpLayout() {

        reference.setText(properties.getProperty("ReferenceFile"));
        input.setText(properties.getProperty("InputFile"));
        if (!input.getText().isEmpty())
            setOutputFiles();
        cores.setText(String.valueOf(Runtime.getRuntime().availableProcessors()));
        cores.setTextFormatter(coresFormatter);
        multiplicity.setText(String.valueOf(MyParameters.MULTIPLICITY));
        multiplicity.setTextFormatter(multiplicityFormatter);

        ObservableList<String> typeList = FXCollections.observableArrayList(READ_TYPES);
        typeBox.setItems(typeList);
        typeBox.getSelectionModel().selectFirst();

        matrices.setItems(matrixList);
        matrices.getSelectionModel().selectFirst();
        adaptPenalties();
        frameshift.setText(String.valueOf(MyParameters.FRAMESHIFT_PENALTY));
        frameshift.setTextFormatter(frameshiftFormatter);

        minCoverage.setValue(MyParameters.MIN_COVERAGE);
        minCovLabel.setText(((int) minCoverage.getValue()) + "%");
        eValue.setText(String.valueOf(MyParameters.MAX_EVALUE));
        eValue.setTextFormatter(evalueFormatter);

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

    private void chooseFile(StringProperty s, String title, String key) {
        FileChooser chooser = new FileChooser();
        chooser.setTitle(title);
        File dir = chooser.showOpenDialog(stage);
        if (dir != null) {
            s.set(dir.getAbsolutePath());
            setProperties(key, dir.getAbsolutePath());
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
                        String fileName = new File(f).getName().replaceAll(ending, ".daa");
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

    private File checkFile(TextField text) {
        if (text.getText() == null || text.getText().isEmpty()) {
            info.setText("WARNING: file not specified!");
            return null;
        }
        File file = new File(text.getText());
        if (file.exists())
            return file;
        info.setText("WARNING: " + text.getText() + " is not a file!");
        return null;
    }

    private File checkFolder(TextField text) {
        if (text.getText() == null || text.getText().isEmpty()) {
            info.setText("WARNING: file not specified!");
            return null;
        }
        File file = new File(text.getText());
        if (file.getParentFile() != null && file.getParentFile().exists())
            return file;
        info.setText("WARNING: " + text.getText() + " is not a folder!");
        return null;
    }

    public void setPresenter(Presenter presenter) {
        this.presenter = presenter;
    }

    public void show() {
        stage.show();
    }
}
