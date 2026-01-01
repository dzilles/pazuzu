import pytest
import warp as wp

def get_available_devices():
    devices = ["cpu"]
    # wp.init() is safe to call multiple times
    wp.init()
    if wp.is_cuda_available():
        devices.append("cuda")
    return devices

@pytest.fixture(params=get_available_devices())
def device(request):
    return request.param
