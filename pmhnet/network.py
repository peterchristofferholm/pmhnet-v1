import numpy as np
import torch
import torch.nn as nn
import torch.nn.functional as F


class SurvProbability(nn.Module):

    def __init__(self, in_features, out_features):
        super().__init__()
        self.out = nn.Linear(in_features, out_features)

    def forward(self, x):
        hazard = torch.sigmoid(self.out(x))
        probs = torch.stack([1 - hazard, hazard], dim=1)
        return probs


class NegLogLik(nn.Module):

    def __init__(self):
        super().__init__()

    def forward(self, predicted, observed):
        eps = 1e-7
        batch_size = observed.shape[0]
        negloglik = torch.sum(-torch.log(predicted[observed] + eps))
        return negloglik / batch_size


class HazToSurv(nn.Module):
    """
    Converts a tensor of conditional hazards to a tensor with survival
    probabilities. A prediction for an exact time `xout` is approximated using
    linear interpolation.

    breaks: with 0.0 included at index 0
    """

    def __init__(self, breaks, xout):

        super().__init__()
        self.i = np.argmax(breaks > xout) - 1
        self.x1 = breaks[self.i]
        self.x2 = breaks[self.i+1]
        self.xout = xout

    def forward(self, ys):

        ys = torch.cumprod(ys[:, 0], dim=1)
        y1, y2 = ys[:, self.i], ys[:, self.i+1]
        surv = ((y2-y1)/(self.x2-self.x1)) * (self.xout-self.x1) + y1
        return torch.unsqueeze(surv, 1)
