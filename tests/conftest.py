"""Shared test setup: the suite runs without network access (issue #19).

Every test gets an autouse guard that refuses any socket connection or name lookup
outside the machine. A test that reaches ENA, NCBI or AWS by accident fails with
NetworkBlocked instead of passing on a good day and failing on a bad one. Tests that
need a server use a fake (monkeypatched requests.get, routes._run, ...) and never
open a socket, so the guard does not touch them.

Loopback and Unix sockets stay allowed: they do not leave the machine.

The guard covers this Python process only. A test that starts a subprocess is not
covered; those tests run local scripts against committed files.

Set POPSEQ_ALLOW_NETWORK=1 to switch the guard off for a deliberate live test.
"""
from __future__ import annotations

import os
import socket

import pytest

_LOCAL_NAMES = {"localhost", "127.0.0.1", "::1", "", None}


class NetworkBlocked(RuntimeError):
    """A test tried to reach the network."""


def _is_local(address) -> bool:
    if isinstance(address, (str, bytes)):  # AF_UNIX path
        return True
    host = address[0] if isinstance(address, tuple) and address else None
    if isinstance(host, bytes):
        host = host.decode()
    return host in _LOCAL_NAMES or str(host).startswith("127.")


@pytest.fixture(autouse=True)
def _no_network(monkeypatch):
    if os.environ.get("POPSEQ_ALLOW_NETWORK") == "1":
        yield
        return
    real_connect = socket.socket.connect
    real_connect_ex = socket.socket.connect_ex
    real_getaddrinfo = socket.getaddrinfo

    def connect(self, address):
        if not _is_local(address):
            raise NetworkBlocked(f"test tried to connect to {address!r}")
        return real_connect(self, address)

    def connect_ex(self, address):
        if not _is_local(address):
            raise NetworkBlocked(f"test tried to connect to {address!r}")
        return real_connect_ex(self, address)

    def getaddrinfo(host, *args, **kwargs):
        if not _is_local((host,)):
            raise NetworkBlocked(f"test tried to resolve {host!r}")
        return real_getaddrinfo(host, *args, **kwargs)

    monkeypatch.setattr(socket.socket, "connect", connect)
    monkeypatch.setattr(socket.socket, "connect_ex", connect_ex)
    monkeypatch.setattr(socket, "getaddrinfo", getaddrinfo)
    yield
