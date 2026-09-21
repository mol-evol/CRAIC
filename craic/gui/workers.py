"""Tiny async helper so long alignment/analysis work never freezes the UI."""

from __future__ import annotations

import threading
import traceback
from typing import Callable, Set

from PySide6.QtCore import QObject, QRunnable, QThreadPool, Signal, Slot


class _Signals(QObject):
    done = Signal(object)
    error = Signal(str)


class _Task(QRunnable):
    def __init__(self, fn: Callable):
        super().__init__()
        self.fn = fn
        self.signals = _Signals()

    @Slot()
    def run(self):
        try:
            result = self.fn()
        except Exception:
            self._emit(self.signals.error, traceback.format_exc())
            return
        self._emit(self.signals.done, result)

    @staticmethod
    def _emit(signal, payload) -> None:
        # If the signal object was already torn down, there is nothing to
        # deliver — swallow it rather than crash the worker thread.
        try:
            signal.emit(payload)
        except RuntimeError:
            pass


# Strong references to in-flight tasks. Without this, Python can garbage-collect
# the QRunnable and its signal object while the worker thread is still running —
# which silently drops the result (Linux) or crashes with a segfault (macOS, e.g.
# after clicking or double-clicking a residue, which fires an async probe).
_PENDING: Set[QRunnable] = set()


def run_async(fn: Callable, on_done: Callable, on_error: Callable = None) -> None:
    task = _Task(fn)
    _PENDING.add(task)

    def _finish(callback, arg):
        _PENDING.discard(task)        # release once delivered, on the GUI thread
        if callback is not None:
            callback(arg)

    task.signals.done.connect(lambda r: _finish(on_done, r))
    task.signals.error.connect(lambda m: _finish(on_error, m))
    QThreadPool.globalInstance().start(task)


# --------------------------------------------------------------------------- #
# Cancellable variant with progress reporting
# --------------------------------------------------------------------------- #

class CancelToken:
    """Thread-safe stop flag shared between the GUI and a worker."""

    def __init__(self):
        self._event = threading.Event()

    def cancel(self) -> None:
        self._event.set()

    def cancelled(self) -> bool:
        return self._event.is_set()


class _ProgressSignals(QObject):
    done = Signal(object)
    error = Signal(str)
    progress = Signal(int, int, str)        # done, total, label


class _ProgressTask(QRunnable):
    def __init__(self, fn: Callable, token: "CancelToken"):
        super().__init__()
        self.fn = fn
        self.token = token
        self.signals = _ProgressSignals()

    @Slot()
    def run(self):
        def report(done, total, label=""):
            try:
                self.signals.progress.emit(int(done), int(total), str(label))
            except RuntimeError:
                pass
        try:
            result = self.fn(report, self.token.cancelled)
        except Exception:
            try:
                self.signals.error.emit(traceback.format_exc())
            except RuntimeError:
                pass
            return
        try:
            self.signals.done.emit(result)
        except RuntimeError:
            pass


def run_cancellable(fn: Callable, on_done: Callable, on_error: Callable = None,
                    on_progress: Callable = None) -> CancelToken:
    """Run ``fn(report, cancelled)`` off the UI thread.

    ``report(done, total, label)`` drives a progress indicator; ``cancelled()``
    returns True once the caller requests a stop (the work must check it). Returns
    a :class:`CancelToken` whose ``.cancel()`` requests that stop.
    """
    token = CancelToken()
    task = _ProgressTask(fn, token)
    _PENDING.add(task)

    def _finish(callback, arg):
        _PENDING.discard(task)
        if callback is not None:
            callback(arg)

    task.signals.done.connect(lambda r: _finish(on_done, r))
    task.signals.error.connect(lambda m: _finish(on_error, m))
    if on_progress is not None:
        task.signals.progress.connect(lambda d, t, lbl: on_progress(d, t, lbl))
    QThreadPool.globalInstance().start(task)
    return token
