import optuna
from torch.utils.data import DataLoader, ConcatDataset, random_split
from numba import jit

class EarlyStopping:
    """EarlyStopping class for training loop.

    Implementation adapted from [1].

    Usage example:

        es = EarlyStopping(patience=5)
        for epoch in range(num_epochs):
            train_epoch(model, data_loader)
            metric = eval(model, data_loader_dev)
            if es.step(metric):
                break

    [1] https://gist.github.com/stefanonardo/693d96ceb2f531fa05db530f3e21517d
    """

    def __init__(self, patience=100, min_delta=0):
        """
        Arguments:
        ----------
            patience : int
                Number of epochs with no improvement after which training
                should be stopped. Default = 10.
            min_delta : float
                Minimum change in monitored quantity to qualify as an
                improvement. Default = 0.
        """

        self.min_delta = min_delta
        self.patience = patience
        self.best = None  # best metric so far
        self.wait = 0     # epochs without improvement

    def _improvement(self, new, old):

        return new < old - self.min_delta

    def step(self, metric):

        if self.best is None:
            self.best = metric
            return False

        if self._improvement(metric, self.best):
            self.wait = 0
            self.best = metric
        else:
            self.wait += 1

        return self.wait >= self.patience


def split_loaders(dataset, ratio, **params):
    """Splits a dataset using ratio and returns independent dataloaders"""

    length1 = int(len(dataset) * ratio)
    length2 = len(dataset) - length1

    split1, split2 = random_split(dataset, [length1, length2])
    return DataLoader(split1, **params), DataLoader(split2, **params)


def kfold_loaders(dataset, k=5, **params):

    lengths = [len(dataset) // k] * k
    lengths[k-1] += len(dataset) % k

    splits = random_split(dataset, lengths)

    for i in range(k):

        valid_dl = DataLoader(splits[i], **params)

        train_dl = DataLoader(
            ConcatDataset(splits[:i] + splits[i+1:]), **params
        )

        yield train_dl, valid_dl


def reset_weights(model):

    for name, module in model.named_children():
        try:
            module.reset_parameters()
        except AttributeError:
            for subname, submodule in module.named_children():
                submodule.reset_parameters()


@jit(nopython=True)
def concordance_index(time, event, risk):
    """Computes harrels C-index

    Arguments:
        time: 1d follow-up time
        event: 1d indicator if event or censoring
        risk: 1d estimated risk
    """

    assert len({len(time), len(risk), len(event)}) == 1, (
        "All inputs does not have the same length!"
    )

    n_conc, n_disc = 0, 0

    N = len(time)

    for i in range(0, N-1):
        for j in range(i+1, N):

            if event[i] and event[j]:
                if risk[i] > risk[j] and time[i] < time[j]:
                    n_conc += 1
                else:
                    n_disc += 1

            elif event[i]:
                if time[i] < time[j]:
                    if risk[i] > risk[j]:
                        n_conc += 1
                    else:
                        n_disc += 1

            elif event[j]:
                if time[i] > time[j]:
                    if risk[i] < risk[j]:
                        n_conc += 1
                    else:
                        n_disc += 1

    return n_conc / (n_conc + n_disc)


def pred_surv(y_pred, breaks, fut):

    """Predicted survival probability from nnet-survival model

    Arguments:
        y_pred: (n_individuals, n_intervals)
            the conditional probability for surviving each time interval
        breaks: breaks for time-bins starting with 0
        fut: follow-up time point at which predictions are needed

    Returns:
        predicted survival probability for each individual
    """

    y_pred = np.cumprod(y_pred, axis=1)
    n_pers = y_pred.shape[0]

    pred = [np.interp(fut, breaks[1:], y_pred[i, :]) for i in range(n_pers)]

    return np.array(pred)


def study_exists(study_name, storage):

    try:
        _ = optuna.load_study(study_name, storage)
    except KeyError:
        return False

    return True
