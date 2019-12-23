package ella.view.detailsView;

import ella.model.io.MyParameters;
import javafx.beans.property.SimpleStringProperty;
import javafx.collections.FXCollections;
import javafx.collections.ObservableList;
import javafx.geometry.Insets;
import javafx.scene.control.TableColumn;
import javafx.scene.control.TableView;
import javafx.scene.control.cell.PropertyValueFactory;

public class SettingsView extends TableView {

    private ObservableList<Entry> entries;

    private static String[][] alignParam2value = {
            {"--db/-d \t path to the indexed protein reference database file in EDB format", ""},
            {"--in/-i \t path to the query read file in FASTA/FASTQ format", ""},
            {"--out/-o \t path of the reported output file in DAA format", ""},
            {"--threads/-p \t number of parallel threads (default: #cpus)", MyParameters.CPUS + ""},
            {"-F \t penalty for switching frames during blastx alignments", MyParameters.FRAMESHIFT_PENALTY + ""},
            {"-m \t maximum initial matches per query position", MyParameters.MULTIPLICITY + ""},
            {"-e \t maximum e-value for reported alignments", MyParameters.MAX_EVALUE + ""},
            {"-c \t minimum percentage of reference positions covered by reported alignments", MyParameters.MIN_COVERAGE + ""},
            {"--matrix \t scoring matrix", MyParameters.SCORING_MATRIX.getType()},
            {"--gop \t gap open penalty", MyParameters.SCORING_MATRIX.getGapOpen() + ""},
            {"--gep \t gap extension penalty", MyParameters.SCORING_MATRIX.getGapExtend() + ""}
    };

    private static String[][] indexParam2value = {
            {"--in/-i \t path to the input protein reference database file in FASTA format", ""},
            {"--db/-d \t path to the output indexed reference database file in EDB format", ""},
            {"--size/-s \t maximal size of each volume", ""},
            {"--threads/-p \t number of parallel threads (default: #cpus)", ""}
    };

    public void init(String mode, String settings, DetailsView detailsView) {

        setPadding(new Insets(5, 5, 5, 5));

        String[][] params = mode.equals("blastx") ? alignParam2value : indexParam2value;
        String[] split = settings.split("\\s+");
        for (int i = 0; i < split.length; i++) {
            String option = split[i];
            String value = split[++i];
            for (String[] o : params) {
                if (o[0].contains(option))
                    o[1] = value;
            }
        }

        entries = FXCollections.observableArrayList();
        for (String[] o : params) {
            String parameter = o[0].split("\t")[0];
            String description = o[0].split("\t")[1];
            String value = o[1];
            entries.add(new Entry(parameter, value, description));
        }

        TableColumn descriptionColumn = new TableColumn("Description");
        TableColumn parameterColumn = new TableColumn("Parameter");
        TableColumn valueColumn = new TableColumn("Value");
        getColumns().setAll(parameterColumn, valueColumn, descriptionColumn);

        parameterColumn.prefWidthProperty().bind(detailsView.widthProperty().multiply(1. / 7.));
        valueColumn.prefWidthProperty().bind(detailsView.widthProperty().multiply(2./ 7.));
        descriptionColumn.prefWidthProperty().bind(detailsView.widthProperty().multiply(4. / 7.));

        parameterColumn.setCellValueFactory(new PropertyValueFactory<Entry, String>("parameter"));
        valueColumn.setCellValueFactory(new PropertyValueFactory<Entry, String>("value"));
        descriptionColumn.setCellValueFactory(new PropertyValueFactory<Entry, String>("description"));

        setItems(entries);

    }

    public class Entry {

        private SimpleStringProperty parameter;
        private SimpleStringProperty value;
        private SimpleStringProperty description;

        public Entry(String parameter, String value, String description) {
            this.parameter = new SimpleStringProperty(parameter);
            this.value = new SimpleStringProperty(value);
            this.description = new SimpleStringProperty(description);
        }

        public String getParameter() {
            return parameter.get();
        }

        public SimpleStringProperty parameterProperty() {
            return parameter;
        }

        public String getValue() {
            return value.get();
        }

        public SimpleStringProperty valueProperty() {
            return value;
        }

        public String getDescription() {
            return description.get();
        }

        public SimpleStringProperty descriptionProperty() {
            return description;
        }
    }

}
