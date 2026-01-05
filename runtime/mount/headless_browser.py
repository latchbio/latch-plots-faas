"""Minimal Playwright wrapper for headless notebook browser."""

from playwright.async_api import async_playwright, Browser, Page, Playwright


class HeadlessBrowser:
    """Manages a headless Chromium browser for running notebook frontend."""

    def __init__(self) -> None:
        self.playwright: Playwright | None = None
        self.browser: Browser | None = None
        self.page: Page | None = None

    async def start(self, notebook_url: str) -> Page:
        """Launch browser and navigate to notebook.

        Args:
            notebook_url: Full URL to the notebook page with auth params.

        Returns:
            The Page object for the loaded notebook.
        """
        self.playwright = await async_playwright().start()
        self.browser = await self.playwright.chromium.launch(headless=True)
        self.page = await self.browser.new_page(viewport={"width": 1280, "height": 800})
        await self.page.goto(notebook_url, wait_until="networkidle")
        return self.page

    async def wait_for_agent_ready(self, timeout_ms: int = 30000) -> None:
        """Wait for frontend to signal agent is connected.

        Args:
            timeout_ms: Maximum time to wait in milliseconds.

        Raises:
            TimeoutError: If agent doesn't become ready within timeout.
        """
        if self.page is None:
            raise RuntimeError("Browser not started. Call start() first.")
        await self.page.wait_for_selector("[data-agent-ready='true']", timeout=timeout_ms)

    async def stop(self) -> None:
        """Close browser and clean up resources."""
        if self.browser:
            await self.browser.close()
            self.browser = None
        if self.playwright:
            await self.playwright.stop()
            self.playwright = None
        self.page = None

