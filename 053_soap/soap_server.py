"""
SOAP Server Demo Application
サンプルの計算サービスとユーザー情報サービスを提供するSOAPサーバー
"""

from spyne import Application, rpc, ServiceBase, Integer, Unicode, Float
from spyne.protocol.soap import Soap11
from spyne.server.wsgi import WsgiApplication
from wsgiref.simple_server import make_server


class CalculatorService(ServiceBase):
    """計算サービス"""

    @rpc(Integer, Integer, _returns=Integer)
    def add(ctx, a, b):
        """2つの整数を足し算"""
        print(f"add({a}, {b}) called")
        return a + b

    @rpc(Integer, Integer, _returns=Integer)
    def subtract(ctx, a, b):
        """2つの整数を引き算"""
        print(f"subtract({a}, {b}) called")
        return a - b

    @rpc(Integer, Integer, _returns=Integer)
    def multiply(ctx, a, b):
        """2つの整数を掛け算"""
        print(f"multiply({a}, {b}) called")
        return a * b

    @rpc(Float, Float, _returns=Float)
    def divide(ctx, a, b):
        """2つの数値を割り算"""
        print(f"divide({a}, {b}) called")
        if b == 0:
            raise ValueError("Division by zero is not allowed")
        return a / b


class UserService(ServiceBase):
    """ユーザー情報サービス"""

    # 簡単なユーザーデータベース（デモ用）
    users = {
        1: {"name": "田中太郎", "email": "tanaka@example.com", "age": 30},
        2: {"name": "佐藤花子", "email": "sato@example.com", "age": 25},
        3: {"name": "鈴木一郎", "email": "suzuki@example.com", "age": 35},
    }

    @rpc(Integer, _returns=Unicode)
    def get_user_name(ctx, user_id):
        """ユーザーIDから名前を取得"""
        print(f"get_user_name({user_id}) called")
        user = UserService.users.get(user_id)
        if user:
            return user["name"]
        return "User not found"

    @rpc(Integer, _returns=Unicode)
    def get_user_email(ctx, user_id):
        """ユーザーIDからメールアドレスを取得"""
        print(f"get_user_email({user_id}) called")
        user = UserService.users.get(user_id)
        if user:
            return user["email"]
        return "User not found"

    @rpc(Integer, _returns=Integer)
    def get_user_age(ctx, user_id):
        """ユーザーIDから年齢を取得"""
        print(f"get_user_age({user_id}) called")
        user = UserService.users.get(user_id)
        if user:
            return user["age"]
        return -1


def create_app():
    """SOAPアプリケーションを作成"""
    application = Application(
        [CalculatorService, UserService],
        tns='http://example.com/soap',
        in_protocol=Soap11(validator='lxml'),
        out_protocol=Soap11()
    )

    return WsgiApplication(application)


if __name__ == '__main__':
    # サーバーの設定
    HOST = 'localhost'
    PORT = 8000

    print(f"Starting SOAP server on http://{HOST}:{PORT}")
    print(f"WSDL: http://{HOST}:{PORT}/?wsdl")
    print("Press Ctrl+C to stop the server\n")

    # WSGIサーバーを起動
    wsgi_app = create_app()
    server = make_server(HOST, PORT, wsgi_app)

    try:
        server.serve_forever()
    except KeyboardInterrupt:
        print("\nServer stopped")
