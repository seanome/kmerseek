import logging

# Set up logger
logger = logging.getLogger(__name__)


def setup_logging(debug_mode: bool):
    """Configure logging based on debug mode."""
    log_level = logging.DEBUG if debug_mode else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )
