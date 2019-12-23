package ella.model.utils;

import ella.presenter.Presenter;

public class RuntimeThread extends Thread {

    private long time = System.currentTimeMillis();
    private boolean isDone = false;
    private Presenter presenter;
    private int taskId;

    public RuntimeThread(int taskId, Presenter presenter) {
        this.taskId = taskId;
        this.presenter = presenter;
    }

    public void run() {
        while (!isDone) {
            try {
                Thread.sleep(1000);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            if (!isDone) {
                String runtime = calculateRuntime((System.currentTimeMillis() - time) / 1000);
                presenter.reportRuntime(taskId, runtime);
            }
        }
    }

    private String calculateRuntime(long timeInMillis) {
        long seconds = (int) ((timeInMillis) % 60);
        long minutes = (int) ((timeInMillis / (60)) % 60);
        long hours = (int) ((timeInMillis / (60 * 60)) % 24);
        String hr = (hours < 10) ? "0" + hours : String.valueOf(hours);
        String mn = (minutes < 10) ? "0" + minutes : String.valueOf(minutes);
        String sc = (seconds < 10) ? "0" + seconds : String.valueOf(seconds);
        return hr + ":" + mn + ":" + sc;
    }

    public void setDone(boolean done) {
        isDone = done;
    }

}

