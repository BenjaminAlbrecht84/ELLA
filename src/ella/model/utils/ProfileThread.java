package ella.model.utils;

import com.sun.management.OperatingSystemMXBean;
import ella.presenter.Presenter;

import java.lang.management.ManagementFactory;
import java.text.SimpleDateFormat;
import java.util.Calendar;

public class ProfileThread extends Thread {

    private OperatingSystemMXBean osBean = ManagementFactory.getPlatformMXBean(
            OperatingSystemMXBean.class);
    private SimpleDateFormat sdf = new SimpleDateFormat("HH:mm:ss");
    private Runtime runtime = Runtime.getRuntime();
    private boolean isDone = false;
    private int taskId;

    private Presenter presenter;

    public ProfileThread(int taskId, Presenter presenter) {
        this.taskId = taskId;
        this.presenter = presenter;
    }

    public void run() {

        while (!isDone) {
            try {
                Thread.sleep(5000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            if (!isDone) {
                String timeStamp = sdf.format(Calendar.getInstance().getTime());
                double cpuLoad = osBean.getProcessCpuLoad() * 100.;
                double allocatedMemory = ((double) runtime.totalMemory() / Math.pow(10, 9));
                double freeMemory = ((double) runtime.freeMemory() / Math.pow(10, 9));
                double usedMemory = allocatedMemory - freeMemory;
                presenter.reportCpuLoad(taskId, timeStamp, cpuLoad, usedMemory, allocatedMemory);
            }
        }
    }

    public void setDone(boolean done) {
        isDone = done;
    }

}
