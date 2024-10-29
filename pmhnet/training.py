from copy import deepcopy
from pmhnet.utils import EarlyStopping
from statistics import mean
import math


class Trainer():

    def __init__(
        self, model, loss_fn, optimizer_fn, train_dl, test_dl,
        max_epochs=math.inf
    ):

        self.model = model
        self.loss_fn = loss_fn
        self.optimizer = optimizer_fn(model.parameters())
        self.train_dl = train_dl
        self.test_dl = test_dl
        self.max_epochs = max_epochs

        self.epoch = 0
        self.train_loss = 0
        self.val_loss = 0

    def __iter__(self):

        return self

    def __next__(self):

        self.model.train()  # network in training mode
        running_loss = 0

        for i, batch in enumerate(self.train_dl):

            self.optimizer.zero_grad()  # zero parameter gradients
            labels, features = batch

            predicted = self.model(features)
            loss = self.loss_fn(predicted, labels)
            running_loss += loss.item()
            loss.backward()
            self.optimizer.step()

        self.train_loss = running_loss/(i+1)

        self.model.eval()  # network in validation mode
        running_loss = 0.0

        for i, batch in enumerate(self.test_dl):
            labels, features = batch
            predicted = self.model(features)
            loss = self.loss_fn(predicted, labels)
            running_loss += loss.item() * labels.shape[0]

        self.val_loss = running_loss / len(self.test_dl.sampler)
        self.epoch += 1

        if self.epoch > self.max_epochs:
            raise StopIteration

        return self


class CrossValidation():

    def __init__(self, model, kf_splits, loss_fn, optimizer_fn, **kwargs):

        patience = kwargs.get("patience", math.inf)
        self.max_epochs = kwargs.get("max_epochs", math.inf)

        trainers = []
        for train_dl, test_dl in kf_splits:

            model_copy = deepcopy(model)  # different parallel instances
            trainer = Trainer(
                model_copy, loss_fn, optimizer_fn, train_dl, test_dl
            )
            trainers.append(trainer)

        self.trainers = trainers
        self.early_stopping = EarlyStopping(patience=patience)
        self.epoch = 0
        self.loss = 0

    def __iter__(self):

        return self

    def __next__(self):

        losses = [next(trainer).val_loss for trainer in self.trainers]

        self.loss = mean(losses)
        self.epoch += 1

        if self.early_stopping.step(self.loss) or self.epoch > self.max_epochs:
            raise StopIteration

        return self
