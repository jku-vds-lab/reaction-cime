import os
import resource

import uvicorn

dir = os.path.dirname(os.path.realpath(__file__))


def limit_memory(maxsize):
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    resource.setrlimit(resource.RLIMIT_AS, (maxsize, hard))


if __name__ == "__main__":
    # limit_memory(8 * 1024 * 1024 * 1024)  # MB * KiB * B
    uvicorn.run(
        os.path.basename(dir) + ".dev_app:app",
        host="0.0.0.0",
        port=9000,
        reload=True,
        reload_dirs=[dir],
    )
