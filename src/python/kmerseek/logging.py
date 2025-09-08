import logging
import sys

# Set up logger
logger = logging.getLogger(__name__)


def setup_logging(debug_mode: bool):
    """Configure logging based on debug mode."""
    log_level = logging.DEBUG if debug_mode else logging.INFO

    # Clear existing handlers to ensure our configuration takes effect
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
        stream=sys.stderr,  # Send logging to stderr
        force=True,  # Force reconfiguration
    )
